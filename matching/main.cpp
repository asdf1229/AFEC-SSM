#include "matching/run_matching.h"

int main(int argc, char *argv[])
{
    return ssm::run_algorithm_main(argc, argv, ssm::create_algorithm_definition());
}
