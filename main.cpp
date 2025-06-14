#include <filesystem>
#include "src/BenchmarkDriver.hpp"
#include "src/SimdEverywhere.hpp"
#include "src/algos/Comet.hpp"
#include "src/algos/Skipping.hpp"
#include "src/algos/StartsWith.hpp"
#include "src/algos/StdFind.hpp"
#include "src/algos/Memmem.hpp"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
   if (argc != 3) {
      std::cerr << "Usage: " << argv[0] << " <column:file> <like-pattern:str>" << std::endl;
      exit(-1);
   }

   BenchmarkDriver driver;

   // std::find.
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());
   driver.AddEngine(std::make_unique<StdFindEngineFactory>());

   // std::memmem.
   // driver.AddEngine(std::make_unique<MemmemEngineFactory>());
   // driver.AddEngine(std::make_unique<MemmemEngineFactory>());
   // driver.AddEngine(std::make_unique<MemmemEngineFactory>());

   // std::starts_with.
   driver.AddEngine(std::make_unique<StartsWithEngineFactory>());
   driver.AddEngine(std::make_unique<SkippingEngineFactory>());
   driver.AddEngine(std::make_unique<SkippingEngineFactory>());
   driver.AddEngine(std::make_unique<SkippingEngineFactory>());

   // Comet.
   driver.AddEngine(std::make_unique<CometEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());
   driver.AddEngine(std::make_unique<CometEngineFactory>());

   auto file_path = argv[1];
   auto pattern = argv[2];

   // Check if the file exists
   if (!std::filesystem::exists(file_path)) {
      std::cerr << "Error: File '" << file_path << "' does not exist." << std::endl;
      return 1;
   }

   std::cout << "Running: " << pattern << " on " << file_path << std::endl;
   std::cout << "--------" << std::endl;
   driver.LoadBlocks(file_path);
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
