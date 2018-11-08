#pragma once



namespace kahypar {


class Exchanger {
  

 public: 
  Exchanger(MPI_Comm communicator) : 
    _current_best_fitness(std::numeric_limits<int>::max()),
    _maximum_allowed_pushes(),
    _number_of_pushes(),
    _rank(),
    _individual_already_sent_to(),
    _partition_buffer(),
    _request_list(),
    _m_communicator(communicator) {
      int rank;
      MPI_Comm_rank( _m_communicator, &rank);
      _rank = rank;
      int comm_size;
      MPI_Comm_size( _m_communicator, &comm_size);
      if(comm_size > 2) {_maximum_allowed_pushes = ceil(log2(comm_size));}
      else              {_maximum_allowed_pushes = 1;}
      MPI_Datatype MPI_Partition;
      MPI_Type_contiguous(hg.initialNumNodes(), MPI_UINT, &partition);
      MPI_Type_commit(&partition);
      
      _individual_already_sent_to = std::vector<bool>(comm_size);
    }
  
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

      MPI_Barrier( _m_communicator );
      for( unsigned i = 0; i < _request_list.size(); i++) {
              MPI_Cancel( _request_list[i] );
      }
        
      for( unsigned i = 0; i < _request_list.size(); i++) {
        MPI_Status st;
        MPI_Wait( _request_list[i], & st );
        delete[] _partition_buffer[i];
        delete   _request_list[i];
      }
    }          
  
  inline void sendBestIndividual(const Context& context, const Hypergraph& hg, const Population& population) {

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
      

      
      DBG << "Rank " << _rank << " sending to " << new_target << "...";
      std::vector<PartitionID> partition_vector = population.individualAt(population.best()).partition();
      int* partition_map = new int[partition_vector.size()];

      for(int i = 0; i < partition_vector.size(); ++i) {
        partition_map[i] = partition_vector[i];
        std::cout << partition_vector[i] << " ";
      }
      
      MPI_Request* request = new MPI_Request();
      MPI_Isend(partition_map, hg.initialNumNodes(), MPI_INT, new_target, new_target, _m_communicator, request);
      _number_of_pushes++;
      _individual_already_sent_to[new_target] = true;
      _request_list.push_back(request); 
      //_partition_buffer.push_back(partition_vector.data());
      _partition_buffer.push_back(partition_map);
      
      DBG << "Rank " << _rank << " sending to " << new_target << "success";
    }
    for(unsigned i = 0; i < _request_list.size(); ++i) {
      int finished = 0;
      MPI_Status status;
      MPI_Test(_request_list[i], &finished, &status);
      if(finished) {
        std::swap(_request_list[i], _request_list[_request_list.size()-1]);
        std::swap(_partition_buffer[i], _partition_buffer[_request_list.size()-1]);

        delete[] _partition_buffer[_partition_buffer.size() - 1];
        delete   _request_list[_request_list.size() - 1];

        _partition_buffer.pop_back();
        _request_list.pop_back();
      }
    }
  }
  
  inline void receiveIndividual(const Context& context, Hypergraph& hg, Population& population) {
    int flag; 
    MPI_Status st;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, _m_communicator, &flag, &st);
    
    while(flag){
      DBG << "Rank " << _rank << "recieving...";

      int* partition_vector_pointer = new int[hg.initialNumNodes()];

      
      DBG << hg.initialNumNodes() << _rank;
    
      MPI_Status rst;
      MPI_Recv(partition_vector_pointer, hg.initialNumNodes(), MPI_INT, st.MPI_SOURCE, _rank, _m_communicator, &rst); 
      DBG << "until here breakpoint rank " << _rank;
      std::vector<PartitionID> result_individual_partition(partition_vector_pointer, partition_vector_pointer + hg.initialNumNodes());
      DBG << "Breakpoint A" << "Rank " << _rank;
      for(int i = 0; i < result_individual_partition.size(); ++i) {
        std::cout << result_individual_partition[i] << " ";
      }
      DBG << "Breakpoint B" << "Rank " << _rank;
      hg.reset();
      hg.setPartition(result_individual_partition);
      size_t insertion_value = population.insert(Individual(hg, context), context);
      int recieved_fitness = population.individualAt(insertion_value).fitness();
      DBG << "Rank " << _rank << "recieved Individual from" << st.MPI_SOURCE << "with fitness" << recieved_fitness;
      
      
      if(insertion_value == std::numeric_limits<int>::max()) {
        DBG << "INSERTION DISCARDED";
      }
      
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
      DBG << "Rank " << _rank << "recieving Success";
    }
  }
  inline void exchangeInitialPopulations(Population& population, const Context& context, Hypergraph& hg, const int& amount_of_individuals_already_generated) {
    int number_of_exchanges_required = context.evolutionary.population_size - amount_of_individuals_already_generated - 1;
    if(number_of_exchanges_required < 0) {
      number_of_exchanges_required = 0;
    }
    for(int i = 0; i < number_of_exchanges_required; ++i) {
      exchangeIndividuals(population, context, hg);
    }
    
  }
  inline void exchangeIndividuals(Population& population, const Context& context, Hypergraph& hg) {
    
    
    int amount_of_mpi_processes;
    MPI_Comm_size( _m_communicator, &amount_of_mpi_processes);
    std::vector<int> somehow_permutated(amount_of_mpi_processes);
    std::iota (std::begin(somehow_permutated), std::end(somehow_permutated), 0);
    Randomize::instance().shuffleVector(somehow_permutated, somehow_permutated.size());
    int sending_to = somehow_permutated[_rank];
    int recieving_from = 0;
    for(unsigned i = 0; i < somehow_permutated.size(); ++i) {
      if (somehow_permutated[i] == _rank) {
        recieving_from = i;
        break;
      }
    }
    
    
    std::vector<PartitionID> outgoing_partition = population.individualAt(population.randomIndividual()).partition();
    int* recieved_partition_pointer = new int[outgoing_partition.size()];
    MPI_Status st;
    MPI_Sendrecv( outgoing_partition.data(), outgoing_partition.size(), MPI_INT, sending_to, 0, 
                  recieved_partition_pointer, outgoing_partition.size(), MPI_INT, recieving_from, 0, _m_communicator, &st); 
                  
    DBG << "Rank " << _rank << "sending to " << sending_to << "quick_start";
    std::vector<int> recieved_partition_vector(recieved_partition_pointer, recieved_partition_pointer + outgoing_partition.size());
    hg.reset();
    hg.setPartition(recieved_partition_vector);
    population.insert(Individual(hg, context), context);
  }
  
  
  
  
  
  
  
  
  
  inline void collectBestPartition(const Population& population, const Context& context, Hypergraph& hg) {
     std::vector<PartitionID> best_local_partition(population.individualAt(population.best()).partition());
     int best_local_objective = population.individualAt(population.best()).fitness();
     DBG << "Rank " << _rank << " best local objective" << best_local_objective;
     int best_global_objective = 0;
     //Maybee i need a temporary variable, lets debug
     MPI_Allreduce(&best_local_objective, &best_global_objective, 1, MPI_INT, MPI_MIN, _m_communicator);
     
     DBG << "Rank " << _rank << " best local objective" << best_local_objective;
     int broadcast_rank = std::numeric_limits<int>::max();
     int global_broadcaster = 0;
     if(best_local_objective == best_global_objective) {
       broadcast_rank = _rank;
     }
     DBG << "Rank " << _rank << " broadcast_rank" << broadcast_rank;
     MPI_Allreduce(&broadcast_rank, &global_broadcaster, 1, MPI_INT, MPI_MIN, _m_communicator);
     MPI_Bcast(best_local_partition.data(), hg.initialNumNodes(), MPI_INT, global_broadcaster, _m_communicator);
     
     hg.setPartition(best_local_partition);
  }
  
  
  
  
  
  
  
  
  
  
  
  
 private: 
 
  HyperedgeWeight _current_best_fitness;
  int _maximum_allowed_pushes;
  int _number_of_pushes;
  int _rank;
  std::vector<bool> _individual_already_sent_to;
  std::vector<int*> _partition_buffer;
  std::vector<MPI_Request*> _request_list;
  MPI_Comm _m_communicator;
  
  static constexpr bool debug = true;
};



} //namespace kahypar

