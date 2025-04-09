#include "src/BenchmarkDriver.hpp"
#include "src/algos/StdFind.hpp"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
int main(int, char**)
{
   BenchmarkDriver driver;
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());
   driver.LoadBlocks("data/l_comment.csv");
   driver.RunRaw("special");
   driver.RunCompressed("special");

   return 0;
}
// -------------------------------------------------------------------------------------
