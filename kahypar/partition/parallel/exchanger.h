#pragma once


#include "kahypar/partition/parallel/partition_buffer.h"
namespace kahypar {


class Exchanger {
  

 public: 
  Exchanger(MPI_Comm communicator, size_t hypergraph_size) : 
    _current_best_fitness(std::numeric_limits<int>::max()),
    _generator(),
    _maximum_allowed_pushes(),
    _number_of_pushes(),
    _rank(),
    _individual_already_sent_to(),
    _partition_buffer(),
    _MPI_Partition(),
    _m_communicator(communicator) {
    
      MPI_Comm_rank( _m_communicator, &_rank);
      int comm_size;
      MPI_Comm_size( _m_communicator, &comm_size);
      if(comm_size > 2) {
        _maximum_allowed_pushes = ceil(log2(comm_size));
      }
      else {
        _maximum_allowed_pushes = 1;
      }

      _generator.seed(1);
      _individual_already_sent_to = std::vector<bool>(comm_size);
      
      MPI_Type_contiguous(hypergraph_size, MPI_INT, &_MPI_Partition);
      MPI_Type_commit(&_MPI_Partition);
    }
  Exchanger(const Exchanger&) = delete;
  Exchanger& operator= (const Exchanger&) = delete;

  Exchanger(Exchanger&&) = delete;
  Exchanger& operator= (Exchanger&&) = delete;

  ~Exchanger() {
      MPI_Barrier( _m_communicator );
        
      int flag;
      MPI_Status st;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _m_communicator, &flag, &st);
        
      while(flag) {
        int message_length;
        MPI_Get_count(&st, MPI_INT, &message_length);
                 
        int* partition_map = new int[message_length];
        MPI_Status rst;
        MPI_Recv( partition_map, message_length, MPI_INT, st.MPI_SOURCE, _rank, _m_communicator, &rst); 
                
        delete[] partition_map;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _m_communicator, &flag, &st);
      } 
      MPI_Type_free(&_MPI_Partition);
      MPI_Barrier( _m_communicator );
    }          
  inline void clearBuffer() {
    _partition_buffer.releaseBuffer();
  }
  inline void sendBestIndividual(const Population& population) {

    if(population.individualAt(population.best()).fitness() < _current_best_fitness) {
      
      _current_best_fitness = population.individualAt(population.best()).fitness();
      
      
      for(unsigned i = 0; i < _individual_already_sent_to.size(); ++i) {
        _individual_already_sent_to[i] = false;
      }  
      _individual_already_sent_to[_rank] = true;
      _number_of_pushes = 0;
    }
    
    bool something_todo = false;
    for(unsigned i = 0; i < _individual_already_sent_to.size(); ++i) {
      if(!_individual_already_sent_to[i]) {
        something_todo = true;
        break;
      }
    }
    if(_number_of_pushes > _maximum_allowed_pushes) {
      something_todo = false;
    }
    
    if(something_todo) {
      int new_target = _rank;
      
      
      //Determining the target via randomization
      std::vector<int> randomized_targets(_individual_already_sent_to.size());
      std::iota(std::begin(randomized_targets), std::end(randomized_targets), 0);
      Randomize::instance().shuffleVector(randomized_targets, randomized_targets.size());
      
      for(unsigned i = 0; i < _individual_already_sent_to.size(); ++i) { 
        int current_target = randomized_targets[i];
	      if(!_individual_already_sent_to[current_target]) {
          new_target = current_target;
	        break;
        }
      }
      

      
      DBG << "Rank " << _rank << " sending to " << new_target << "..."<< "fitness " << population.individualAt(population.best()).fitness();
      const std::vector<PartitionID>& partition_vector = population.individualAt(population.best()).partition();

      BufferElement ele = _partition_buffer.acquireBuffer(partition_vector);
      DBG << "Rank " << _rank << " buffersize: " << _partition_buffer.size();
      MPI_Isend(ele.partition, 1, _MPI_Partition, new_target, new_target, _m_communicator, ele.request);
      _number_of_pushes++;
      _individual_already_sent_to[new_target] = true;
      //_partition_buffer.acquireBuffer(request, partition_map);
    }
    //_partition_buffer.releaseBuffer();

  }
  
  inline void receiveIndividual(const Context& context, Hypergraph& hg, Population& population) {

    int flag; 
    MPI_Status st;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _m_communicator, &flag, &st);
    
    while(flag) {
    int* partition_vector_pointer = new int[hg.initialNumNodes()];
    MPI_Status rst;
    MPI_Recv(partition_vector_pointer, 1, _MPI_Partition, st.MPI_SOURCE, _rank, _m_communicator, &rst); 
    std::vector<PartitionID> result_individual_partition(partition_vector_pointer, partition_vector_pointer + hg.initialNumNodes());
    delete[] partition_vector_pointer;
    hg.reset();
    hg.setPartition(result_individual_partition);
    size_t insertion_value = population.insert(Individual(hg, context), context);
    if(insertion_value == std::numeric_limits<unsigned>::max()) {
      DBG << "INSERTION DISCARDED";
      return;
    }
    else {
      LOG <<" MPIRank " << context.mpi.rank << ":"  << "Population " << population << "receive individual exchanger.h l.148";
    }
    int recieved_fitness = population.individualAt(insertion_value).fitness();
    DBG << "Rank " << _rank << "recieved Individual from" << st.MPI_SOURCE << "with fitness" << recieved_fitness;
    
    DBG << "Rank " << _rank << " buffersize: " << _partition_buffer.size(); 
    if(recieved_fitness < _current_best_fitness) {
      
      _current_best_fitness = recieved_fitness;
                       
      for(unsigned i = 0; i < _individual_already_sent_to.size(); ++i) {
        _individual_already_sent_to[i] = false;
      }
      _individual_already_sent_to[_rank] = true;
      _number_of_pushes = 0;
    }
    _individual_already_sent_to[st.MPI_SOURCE] = true;
      
     MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _m_communicator, &flag, &st);
     }
  }
  
  inline void exchangeInitialPopulations(Population& population, const Context& context, Hypergraph& hg, const int& amount_of_individuals_already_generated) {
    int number_of_exchanges_required = context.evolutionary.population_size - amount_of_individuals_already_generated;
    if(number_of_exchanges_required < 0) {
      number_of_exchanges_required = 0;
    }
    for(int i = 0; i < number_of_exchanges_required; ++i) {
      exchangeIndividuals(population, context, hg);
      DBG << "Population " << population << "Rank " << _rank;
    }
    
  }
  inline void exchangeIndividuals(Population& population, const Context& context, Hypergraph& hg) {

    
    int amount_of_mpi_processes;
    MPI_Comm_size( _m_communicator, &amount_of_mpi_processes);
    std::vector<int> permutation_of_mpi_process_numbers(amount_of_mpi_processes);
    DBG << "Rank " << _rank << "length " << permutation_of_mpi_process_numbers.size();

    

    /*This randomization is independent from the algorithm and has to use its own seed to
      ensure that all Mpi threads generate the same permutation */  
      
        
    std::iota (std::begin(permutation_of_mpi_process_numbers), std::end(permutation_of_mpi_process_numbers), 0);
    std::shuffle(permutation_of_mpi_process_numbers.begin(),
                 permutation_of_mpi_process_numbers.begin() + permutation_of_mpi_process_numbers.size() , _generator);
    
    int sending_to = permutation_of_mpi_process_numbers[_rank];
    int recieving_from = 0;
    for(unsigned i = 0; i < permutation_of_mpi_process_numbers.size(); ++i) {
      if (permutation_of_mpi_process_numbers[i] == _rank) {
        recieving_from = i;
        break;
      }
    }
       
    std::vector<PartitionID> outgoing_partition = population.individualAt(population.randomIndividual()).partition();
    int* outgoing_partition_map = new int[outgoing_partition.size()];
    for(size_t i = 0; i < outgoing_partition.size(); ++i) {
      outgoing_partition_map[i] = outgoing_partition[i];
    }
    int* recieved_partition_pointer = new int[outgoing_partition.size()];
    DBG << "Rank " << _rank << "sending to " << sending_to << "quick_start";
    MPI_Status st;
    MPI_Sendrecv( outgoing_partition.data(), 1, _MPI_Partition, sending_to, 0, 
                  recieved_partition_pointer, 1, _MPI_Partition, recieving_from, 0, _m_communicator, &st); 
                  

    std::vector<int> recieved_partition_vector(recieved_partition_pointer, recieved_partition_pointer + outgoing_partition.size());
    
    
    hg.reset();
    hg.setPartition(recieved_partition_vector);
    population.insert(Individual(hg, context), context);
    LOG <<" MPIRank " << context.mpi.rank << ":"  << "Population " << population << "exchange individuals exchanger.h l.225";
    delete[] outgoing_partition_map;
    delete[] recieved_partition_pointer;
  }
  
  
  
  
  
  
  
  
  
  inline void collectBestPartition(Population& population, Hypergraph& hg, const Context& context) {
    DBG << "Collect Best Partition" << " Rank " << _rank;
    std::vector<PartitionID> best_local_partition(population.individualAt(population.best()).partition());
    int best_local_objective = population.individualAt(population.best()).fitness();
    int best_global_objective = 0;
    
    MPI_Allreduce(&best_local_objective, &best_global_objective, 1, MPI_INT, MPI_MIN, _m_communicator);
    
    int broadcast_rank = std::numeric_limits<int>::max();
    int global_broadcaster = 0;
    DBG <<"Rank " << _rank << " " << best_local_objective << " " << best_global_objective;
    if(best_local_objective == best_global_objective) {
      broadcast_rank = _rank;
    }
    
    MPI_Allreduce(&broadcast_rank, &global_broadcaster, 1, MPI_INT, MPI_MIN, _m_communicator);
    DBG << "Rank " << _rank << " " << broadcast_rank;
    MPI_Bcast(best_local_partition.data(), hg.initialNumNodes(), MPI_INT, global_broadcaster, _m_communicator);
     
    hg.setPartition(best_local_partition);
    population.insert(Individual(hg, context), context);
  }
  

  
 private: 
 
  
  HyperedgeWeight _current_best_fitness;
  std::mt19937 _generator;
  int _maximum_allowed_pushes;
  int _number_of_pushes;
  int _rank;
  std::vector<bool> _individual_already_sent_to;
  PartitionBuffer _partition_buffer;
  MPI_Datatype _MPI_Partition;
  MPI_Comm _m_communicator;
  
  static constexpr bool debug = true;
};



} //namespace kahypar

