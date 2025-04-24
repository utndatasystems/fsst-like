#include "src/BenchmarkDriver.hpp"
#include "src/algos/StartsWith.hpp"
#include "src/algos/StdFind.hpp"
#include "src/algos/Comet.hpp"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
int main(int, char**)
{
   BenchmarkDriver driver;
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());
   driver.AddEngine(std::make_unique<StartsWithEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());

   std::cout << "contains" << std::endl;
   std::cout << "--------" << std::endl;
   driver.LoadBlocks("data/l_comment.csv");
   driver.Run("%special%");

   std::cout << "" << std::endl;
   std::cout << "prefix" << std::endl;
   std::cout << "------" << std::endl;
   driver.LoadBlocks("data/p_type.csv");
   driver.Run("MEDIUM POLISHED%"); // not like
   driver.Run("MEDIUM POLISHED%"); // not like
   driver.Run("MEDIUM POLISHED%"); // not like
   // driver.Run("PROMO%");

   return 0;
}
// -------------------------------------------------------------------------------------
