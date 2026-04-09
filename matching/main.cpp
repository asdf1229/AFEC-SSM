#include "matching/run_matching.h"

int main(int argc, char *argv[])
{
    return ssm_ged::run_algorithm_main(argc, argv, ssm_ged::create_algorithm_definition());
}
