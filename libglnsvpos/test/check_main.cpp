//#include <gtest/gtest.h>
//#include <stdlib.h>

//#include "check_suites.h"

//int main(void)
//{
//  int number_failed;
//  Suite *s = position_suite();
//  SRunner *sr = srunner_create(s);
//  srunner_set_xml(sr, "test_results.xml");

//  srunner_set_fork_status(sr, CK_NOFORK);
//  srunner_run_all(sr, CK_NORMAL);
//  number_failed = srunner_ntests_failed(sr);
//  srunner_free(sr);
//  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
//}
