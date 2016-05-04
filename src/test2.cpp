/* -------------------------------------------------------------------------- *
 *
 * IJRR 2013 tests: compare the cost of the HQP with the QP cascade of Kanoun
 * and DeLassa.
 *
 * -------------------------------------------------------------------------- */

#include "soth/HCOD.hpp"




#include <boost/assign/std/vector.hpp> // for 'operator+=()'
using namespace boost::assign; // bring 'operator+=()' into scope


#include <iostream>
#include <sstream>
#include <kdl/expressiontree.hpp>
#include <kdl/conversions.hpp>
#include <kdl/frames_io.hpp>
using namespace soth;

using namespace KDL;
using std::endl;
using std::cout;


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
	//compute jacobian
	//robot in expression
	Expression<Vector>::Ptr L      = Constant(Vector(0,0,0.2));
	Expression<Vector>::Ptr L0=  Constant(Vector(0,0,0.0));
	Expression<Rotation>::Ptr R0=  Constant(Rotation(KDL::Rotation::Identity()));

	// calculate the forward kinematics
	double arm_length = 0.45;
   double forearm_length = 0.45;
   double tool_length = 0.2;
   // get the expression frames
   Expression < Frame >::Ptr w_T_1 = KDL::Constant( KDL::Frame::DH_Craig1989( 0, 0, 0, 0 ) )
									 * KDL::frame( rot_z( KDL::input( 0 ) ) );
   Expression < Frame >::Ptr w_T_2 = w_T_1 * KDL::Constant( KDL::Frame::DH_Craig1989( 0, -M_PI / 2., 0, 0 ) )
									 * KDL::frame( rot_z( KDL::input( 1 ) ) );
   Expression < Frame >::Ptr w_T_3 = w_T_2 * KDL::Constant( KDL::Frame::DH_Craig1989( arm_length, 0, 0, -M_PI / 2. ) )
									 * KDL::frame( rot_z( KDL::input( 2 ) ) );
   Expression < Frame >::Ptr w_T_4= w_T_3 * KDL::Constant( KDL::Frame::DH_Craig1989( 0, -M_PI / 2., forearm_length, 0 ) )
									* KDL::frame( rot_z( KDL::input( 3 ) ) );
   Expression < Frame >::Ptr w_T_5 = w_T_4 * KDL::Constant( KDL::Frame::DH_Craig1989( 0, M_PI / 2., 0, M_PI / 2. ) )
									 * KDL::frame( rot_z( KDL::input( 4 ) ) );
   Expression < Frame >::Ptr w_T_ee = w_T_5 * KDL::Constant( KDL::Frame::DH_Craig1989( 0, M_PI / 2., 0, 0 ) )
									  * KDL::frame( rot_z( KDL::input( 5 ) ) )
									  * KDL::Constant( KDL::Frame(KDL::Vector(0,0,tool_length)));

	//Expression<Rotation>::Ptr w_R_ee=rotation(w_T_ee);
	Expression<Vector>::Ptr w_P_ee=(origin( (w_T_ee)));
	Expression<double>::Ptr w_Px_4=coord_x(origin( (w_T_4)));
	unsigned int NB_STAGE,NC;
	std::vector<unsigned int> NR;
	std::vector<Eigen::MatrixXd> J;
	std::vector<soth::VectorBound> b;



	NB_STAGE = 2;

	NR  += 3,1;
	NC = 6;

	std::vector<int> joint_indexes(NC);
	joint_indexes+=0,1,2,3,4,5;


	std::vector<double> q(NC);
	q+=1.1,1.2,1.3,1.4,1.5,1.6;


	std::cout << "nVar \t= " << NC << std::endl;
	std::cout << "nLevels \t= " << NB_STAGE << std::endl;
	std::cout << "LevelDim \t= [ ";
	for(int i=0;i<(int)NB_STAGE;++i) std::cout << NR[i] << " ";
	std::cout << "]" << std::endl;

	//generateDeficientDataSet(J,b,NB_STAGE,RANKFREE,RANKLINKED,NR,NC);

	HCOD hsolver(NC,NB_STAGE); VectorXd solution(NC);
	w_T_ee->setInputValues(joint_indexes,q);

	b.resize(NB_STAGE);
	J.resize(NB_STAGE);
	cout<<" w_P_ee: "<<w_P_ee->value()<<endl;//always compute value before jacobian!!!
	cout<<" w_Px_4: "<<w_Px_4->value()<<endl;

	for( unsigned int i=0;i<NB_STAGE;++i )
	{
		J[i].resize(NR[i],NC);
		b[i].resize(NR[i]);
		switch (i) {
		case 0:

			for( int j=0;j<b[i].size();++j ){
				b[i][j] = Bound( (j+1)/10.0,(j+2)/10.0 );//explicit cast to double!
			}
			for(unsigned int iii=0;iii<NC;iii++)
			{
				KDL::Vector tw=w_P_ee->derivative(iii);
				J[i].col(iii)=toEigen(tw);
			}
			break;
		default:
			b[i].resize(NR[i]);
			for( int j=0;j<b[i].size();++j ){
				b[i][j] = Bound( (j+1)/10.0,(j+2)/10.0 );//explicit cast to double!
			}
			for(unsigned int iii=0;iii<NC;iii++)
			{
				double d=w_Px_4->derivative(iii);
				J[i](iii)=d;
			}
			break;

		}
		cout<<"bounds "<<i<<"\n"<<b[i]<<endl;
		cout<<"Jacobian "<<i<<"\n"<<J[i]<<endl;
		hsolver.pushBackStage(J[i],b[i]);
		//hsolver.setNameByOrder("level");
	}
	//hsolver.setDamping(0.0);
	hsolver.setInitialActiveSet();


	for(int i=0;i<1;++i)
	{
		hsolver.activeSearch(solution);
		hsolver.show(cout);
		cout<<"Joint velocity\n"<<solution.transpose()<<endl;
		cout<<"twist 1 \n"<<J[0]*solution<<endl;
		cout<<"twist 2 \n"<<J[1]*solution<<endl;
	}

}

