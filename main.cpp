#include <filesystem>
#include "src/BenchmarkDriver.hpp"
#include "src/SimdEverywhere.hpp"
#include "src/algos/Comet.hpp"
#include "src/algos/Skipping.hpp"
#include "src/algos/StartsWith.hpp"
#include "src/algos/StdFind.hpp"
#include "src/algos/Memmem.hpp"
#include "src/algos/Equals.hpp"
#include "src/algos/OnPairLike.hpp"
#include "src/algos/PassThrough.hpp"
#include "src/algos/TokenLike.hpp"
// -------------------------------------------------------------------------------------
using namespace std;
// -------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
   if (argc < 3 || argc > 4) {
      std::cerr << "Usage: " << argv[0] << " <column:file> <like-pattern:str> [fsst|token|onpair]" << std::endl;
      exit(-1);
   }

   BenchmarkDriver driver;
   CompressionCodec codec = CompressionCodec::Fsst;
   if (argc == 4) {
      std::string codec_arg = argv[3];
      if (codec_arg == "fsst") {
         codec = CompressionCodec::Fsst;
      } else if (codec_arg == "token") {
         codec = CompressionCodec::Tokenizer;
      } else if (codec_arg == "onpair") {
         codec = CompressionCodec::OnPair;
      } else {
         std::cerr << "Unknown codec '" << codec_arg << "'. Expected one of: fsst, token, onpair." << std::endl;
         return 1;
      }
   }

   if (codec == CompressionCodec::OnPair) {
      driver.AddEngine(std::make_unique<OnPairLikeEngineFactory>());
      driver.AddEngine(std::make_unique<OnPairLikeEngineFactory>());
      driver.AddEngine(std::make_unique<OnPairLikeEngineFactory>());
   } else {
      if (codec == CompressionCodec::Tokenizer) {
         driver.AddEngine(std::make_unique<TokenLikeEngineFactory>());
         driver.AddEngine(std::make_unique<TokenLikeEngineFactory>());
         driver.AddEngine(std::make_unique<TokenLikeEngineFactory>());
      }

      // std::find.
      driver.AddEngine(std::make_unique<StdFindEngineFactory>());
      driver.AddEngine(std::make_unique<StdFindEngineFactory>());
      driver.AddEngine(std::make_unique<StdFindEngineFactory>());

      // Equality (separate engine; only active for patterns without wildcards).
      driver.AddEngine(std::make_unique<EqualsEngineFactory>());
      driver.AddEngine(std::make_unique<EqualsEngineFactory>());
      driver.AddEngine(std::make_unique<EqualsEngineFactory>());

      // Lower-bound O(n) row walk baseline.
      driver.AddEngine(std::make_unique<PassThroughEngineFactory>());
      driver.AddEngine(std::make_unique<PassThroughEngineFactory>());
      driver.AddEngine(std::make_unique<PassThroughEngineFactory>());

      // std::memmem.
      // driver.AddEngine(std::make_unique<MemmemEngineFactory>());
      // driver.AddEngine(std::make_unique<MemmemEngineFactory>());
      // driver.AddEngine(std::make_unique<MemmemEngineFactory>());

      if (codec == CompressionCodec::Fsst) {
         // FSST-specific compressed-path experiments.
         // driver.AddEngine(std::make_unique<StartsWithEngineFactory>());
         // driver.AddEngine(std::make_unique<SkippingEngineFactory>());
         // driver.AddEngine(std::make_unique<SkippingEngineFactory>());
         // driver.AddEngine(std::make_unique<SkippingEngineFactory>());

         // Comet.
         driver.AddEngine(std::make_unique<CometEngineFactory>());
         driver.AddEngine(std::make_unique<CometEngineFactory>());
         driver.AddEngine(std::make_unique<CometEngineFactory>());
      }
   }

   auto file_path = argv[1];
   auto pattern = argv[2];

   // Check if the file exists
   if (!std::filesystem::exists(file_path)) {
      std::cerr << "Error: File '" << file_path << "' does not exist." << std::endl;
      return 1;
   }

   std::string codec_name = "fsst";
   if (codec == CompressionCodec::Fsst) {
      codec_name = "fsst";
   } else if (codec == CompressionCodec::Tokenizer) {
      codec_name = "token";
   } else if (codec == CompressionCodec::OnPair) {
      codec_name = "onpair";
   }
   std::cout << "Running: " << pattern << " on " << file_path << " (codec=" << codec_name << ")" << std::endl;
   std::cout << "--------" << std::endl;
   driver.LoadBlocks(file_path, codec);
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
