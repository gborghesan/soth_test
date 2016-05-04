/* -------------------------------------------------------------------------- *
 *
 * IJRR 2013 tests: compare the cost of the HQP with the QP cascade of Kanoun
 * and DeLassa.
 *
 * -------------------------------------------------------------------------- */
#define SOTH_DEBUG
#define SOTH_DEBUG_MODE 45
#include "soth/debug.hpp"
#include "soth/HCOD.hpp"
#include "soth/Random.hpp"
#include "RandomGenerator.hpp"
#include "soth/DestructiveColPivQR.hpp"

#include <boost/assign/std/vector.hpp> // for 'operator+=()'
using namespace boost::assign; // bring 'operator+=()' into scope

#ifndef WIN32
#include <sys/time.h>
#endif // WIN32

#include <iostream>
#include <sstream>

using namespace soth;
using std::endl;
using std::cout;
using std::cerr;

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

void ehqp( HCOD & hsolver )
{
  hsolver.initialize();
  //hsolver.Y.computeExplicitly();
  //hsolver.computeSolution();
  //hsolver.showActiveSet(std::cout);
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


int main ()
{
  unsigned int NB_STAGE,NC;
  std::vector<unsigned int> NR,RANKLINKED,RANKFREE;
  std::vector<Eigen::MatrixXd> J;
  std::vector<soth::VectorBound> b;

  soth::Random::setSeed(704819);
  //  const int size = 100;

  //generateFixedSizeRandomProfile(size,
  //                               1,0.99,0.99,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

	  NB_STAGE = 3;
	  RANKFREE += 3,4,3;
	  RANKLINKED += 0,0,0;
	  NR = RANKFREE;
	  NC = 10;


  std::cout << "nVar \t= " << NC << std::endl;
  std::cout << "nLevels \t= " << NB_STAGE << std::endl;
  std::cout << "LevelDim \t= [ ";
  for(int i=0;i<(int)NB_STAGE;++i) std::cout << NR[i] << " ";
  std::cout << "]" << std::endl;

  generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

  HCOD hsolver(NC,NB_STAGE);
  for( unsigned int i=0;i<NB_STAGE;++i )
	{
	  for( int j=0;j<b[i].size();++j )
	b[i][j] = Bound(rand());
	  cout<<"b"<<b[i]<<endl;
	  cout<<"J"<<J[i]<<endl;
	  hsolver.pushBackStage(J[i],b[i]);
	  hsolver.setNameByOrder("level");
	}
  hsolver.setDamping(0.0);
  hsolver.setInitialActiveSet();


  for(int i=0;i<10;++i) ehqp(hsolver);




}
