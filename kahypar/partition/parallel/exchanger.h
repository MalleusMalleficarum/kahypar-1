#pragma once



namespace kahypar {

class Exchanger {

 private: 
 
  HyperedgeWeight current_best_fitness;
  int maximum_allowed_pushes;
  int number_of_pushes;
  std::vector<bool> individual_already_sent_to;
  std::vector<int*> _partition_buffer;
  std::vector<MPI_Request*> _request_list;
  MPI_Comm m_communicator;
  
  static constexpr bool debug = true;
 public: 
  Exchanger(MPI_Comm communicator) : 
    current_best_fitness(std::numeric_limits<int>::max()),
    m_communicator(communicator) {
  
      int comm_size;
      MPI_Comm_size( m_communicator, &comm_size);
    
      if(comm_size > 2) maximum_allowed_pushes = ceil(log2(comm_size));
      else              maximum_allowed_pushes = 1;
      
      individual_already_sent_to = std::vector<bool>(comm_size);
    }
  
  ~Exchanger() {
      MPI_Barrier( m_communicator );
      int rank;
      MPI_Comm_rank( m_communicator, &rank);
        
      int flag; MPI_Status st;
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
        
      while(flag) {
        int message_length;
        MPI_Get_count(&st, MPI_INT, &message_length);
                 
        int* partition_map = new int[message_length];
        MPI_Status rst;
        MPI_Recv( partition_map, message_length, MPI_INT, st.MPI_SOURCE, rank, m_communicator, &rst); 
                
        delete[] partition_map;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
      } 

      MPI_Barrier( m_communicator );
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
  
  inline void broadcastBestIndividual(const Context& context, const Hypergraph& hg, const Population& population, bool final_collection = false) {

    int rank, size;
    MPI_Comm_rank( m_communicator, &rank);
    MPI_Comm_size( m_communicator, &size);
 

    if(population.individualAt(population.best()).fitness() < current_best_fitness || final_collection) {
      
      current_best_fitness = population.individualAt(population.best()).fitness();
      
      
      for(unsigned i = 0; i < individual_already_sent_to.size(); ++i) {
        individual_already_sent_to[i] = false;
      }  
      individual_already_sent_to[rank] = true;
      number_of_pushes = 0;
    }
    
    bool something_todo = false;
    for(unsigned i = 0; i < individual_already_sent_to.size(); ++i) {
      if(!individual_already_sent_to[i]) {
        something_todo = true;
        break;
      }
    }
    if(number_of_pushes > maximum_allowed_pushes) {
      something_todo = false;
    }
    
    if(something_todo) {
      int new_target = rank;
      
      
      //Determining the target via randomization
      std::vector<int> randomized_targets(individual_already_sent_to.size());
      std::iota(std::begin(randomized_targets), std::end(randomized_targets), 0);
      Randomize::instance().shuffleVector(randomized_targets, randomized_targets.size());
      for(unsigned i = 0; i < individual_already_sent_to.size(); ++i) { 
        int current_target = randomized_targets[i];
	      if(!individual_already_sent_to[current_target]) {
          new_target = current_target;
	        break;
        }
      }
      
      if(final_collection) {
        new_target = 0;
      }
      
      DBG << "Rank " << rank << " sending to " << new_target << "...";
      std::vector<PartitionID> partition_vector = population.individualAt(population.best()).partition();
      int* partition_map = new int[partition_vector.size()];
      for(int i = 0; i < partition_vector.size(); ++i) {
        partition_map[i] = partition_vector[i];
      }

     
      MPI_Request* request = new MPI_Request();
      MPI_Isend(partition_map, hg.initialNumNodes(), MPI_INT, new_target, new_target, m_communicator, request);
      number_of_pushes++;
      individual_already_sent_to[new_target] = true;
      _request_list.push_back(request); 
      _partition_buffer.push_back(partition_map);
      
      DBG << "Rank " << rank << " sending to " << new_target << "success";
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
    int rank;
    MPI_Comm_rank( m_communicator, &rank);
    int flag; 
    MPI_Status st;
    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
    
    while(flag){
      DBG << "Rank " << rank << "recieving...";

      int* partition_map = new int[hg.initialNumNodes()];

      
    
    
      MPI_Status rst;
      MPI_Recv(partition_map, hg.initialNumNodes(), MPI_INT, st.MPI_SOURCE, rank, m_communicator, &rst); 
      std::vector<PartitionID> result_individual_partition;
      for(int i = 0; i < hg.initialNumNodes(); ++i) {

        result_individual_partition.push_back(partition_map[i]);
      }

      hg.reset();
      hg.setPartition(result_individual_partition);
      Individual recieved(hg, context);
      DBG << "Rank " << rank << "recieved Individual from" << st.MPI_SOURCE << "with fitness" << recieved.fitness();
      size_t insertion_value = population.insert(Individual(hg, context), context);
      
      if(insertion_value == std::numeric_limits<int>::max()) {
        DBG << "INSERTION DISCARDED";
      }
      
      if(recieved.fitness() < current_best_fitness) {
      
        current_best_fitness = recieved.fitness();
                       
        for(unsigned i = 0; i < individual_already_sent_to.size(); ++i) {
          individual_already_sent_to[i] = false;
        }
        individual_already_sent_to[rank] = true;
        number_of_pushes = 0;
      }
      individual_already_sent_to[st.MPI_SOURCE] = true;
      
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_communicator, &flag, &st);
      DBG << "Rank " << rank << "recieving Success";
    }
  }

};



} //namespace kahypar

