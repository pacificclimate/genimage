#define size_t int

#define NULL 0

class NcFile {
public:
  enum FileMode {
    ReadOnly, // file exists, open read-only
    Write,    // file exists, open for writing
    Replace,  // create new file, even if already exists
    New   // create new file, fail if already exists
  };
  enum FileFormat {
    Classic,         // netCDF classic format (i.e. version 1 format)
    Offset64Bits     // netCDF 64-bit offset format
  };
  NcFile( const char * path, FileMode = ReadOnly ,
          size_t *chunksizeptr = NULL, // optional tuning parameters
          size_t initialsize = 0,
          FileFormat = Classic ) {}
};

int main() {
  NcFile* f = new NcFile("foo", NcFile::Write);
}
