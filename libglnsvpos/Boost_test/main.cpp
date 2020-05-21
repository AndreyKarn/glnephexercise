#define BOOST_TEST_MODULE Test
#include <boost/test/included/unit_test.hpp>

//key config " --log_level=test_suite --run_test=+test_net"

namespace utf = boost::unit_test;

BOOST_AUTO_TEST_CASE(test_1)
{
  BOOST_TEST(true);
}

BOOST_AUTO_TEST_CASE(test_net,
  * utf::description("requires network"))
{
  BOOST_TEST(true);
}
