#include "src/BenchmarkDriver.hpp"
#include "src/SimdEverywhere.hpp"
#include "src/algos/Comet.hpp"
// #include "src/algos/Skipping.hpp"
#include "src/algos/StartsWith.hpp"
#include "src/algos/StdFind.hpp"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <like-pattern:str>" << std::endl;
      exit(-1);
   }

   BenchmarkDriver driver;
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());
   // driver.AddEngine(std::make_unique<StartsWithEngineFactory>());
   // driver.AddEngine(std::make_unique<SkippingEngineFactory>());
   // driver.AddEngine(std::make_unique<SkippingEngineFactory>());
   // driver.AddEngine(std::make_unique<SkippingEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());

   auto pattern = argv[1];

   std::cout << "Running: " << pattern << std::endl;
   std::cout << "--------" << std::endl;
   driver.LoadBlocks("data/l_comment.csv");
   driver.Run(pattern);

   // std::cout << "" << std::endl;
   // std::cout << "prefix" << std::endl;
   // std::cout << "------" << std::endl;
   // driver.LoadBlocks("data/p_type.csv");
   // driver.Run("MEDIUM POLISHED%"); // not like
   // driver.Run("MEDIUM POLISHED%"); // not like
   // driver.Run("MEDIUM POLISHED%"); // not like
   // driver.Run("PROMO%");

   return 0;
}
// -------------------------------------------------------------------------------------
