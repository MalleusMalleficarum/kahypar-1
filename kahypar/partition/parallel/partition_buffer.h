#pragma once



namespace kahypar {

class PartitionBuffer {

 
 

  static constexpr bool debug = true;
 public: 
  PartitionBuffer() :
    _partition_buffer(),
    _request_buffer() {}
  
  
  PartitionBuffer(const PartitionBuffer&) = delete;
  PartitionBuffer& operator= (const PartitionBuffer&) = delete;

  PartitionBuffer(PartitionBuffer&&) = delete;
  PartitionBuffer& operator= (PartitionBuffer&&) = delete;
  
  ~PartitionBuffer() {
    for( unsigned i = 0; i < _request_buffer.size(); i++) {
        MPI_Cancel( _request_buffer[i] );
    }
        
    for( unsigned i = 0; i < _request_buffer.size(); i++) {
      MPI_Status st;
      MPI_Wait( _request_buffer[i], & st );
      delete[] _partition_buffer[i];
      delete   _request_buffer[i];
    }
  }

  void acquireBuffer(MPI_Request* request, int* partition_map) {
    _partition_buffer.push_back(partition_map);
    _request_buffer.push_back(request);  
  }
  
  void releaseBuffer() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    DBG <<"BufferSize "<< _request_buffer.size() << "Rank " << rank;
    for(unsigned i = 0; i < _request_buffer.size(); ++i) {
      int finished = 0;
      MPI_Status status;
      MPI_Test(_request_buffer[i], &finished, &status);
      if(finished) {
        DBG << "Finished " << i << "Rank " << rank;
        std::swap(_request_buffer[i], _request_buffer[_request_buffer.size()-1]);
        std::swap(_partition_buffer[i], _partition_buffer[_request_buffer.size()-1]);

        delete[] _partition_buffer[_partition_buffer.size() - 1];
        delete   _request_buffer[_request_buffer.size() - 1];

        _partition_buffer.pop_back();
        _request_buffer.pop_back();
      }
    }
  }
 private:
 
 std::vector<int*> _partition_buffer;
 std::vector<MPI_Request*> _request_buffer;
};



} //namespace kahypar

