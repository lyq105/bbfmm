// Using cppunit test the mesh io class.
// co
//
#include "cppunit/extensions/HelperMacros.h"

#include "mesh.h"
#include <string>

class MeshTest: public CppUnit::TestFixture
{

	CPPUNIT_TEST_SUITE( MeshTest );
	
	
	CPPUNIT_TEST( test_load_mesh );
	
	
	CPPUNIT_TEST_SUITE_END();


	public:
	MeshTest(){}
	// 初始化类函数
	void setUp();
	// 后处理函数
	void tearDown();
  
	// 测试用例
	int test_load_mesh() 
	{
		Mesh mesh2d,mesh3d;
		mesh2d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat"));
		mesh2d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\2d.dat") + ".plt");

		//mesh2d.load_mesh(std::string("e:\\bem2d_big.dat"));
		//mesh2d.plot_mesh(std::string("e:\\bem2d_big.dat") + ".plt");

		mesh3d.load_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat"));
		mesh3d.plot_mesh(std::string("F:\\LiYiQiang\\Desktop\\Fast Bem\\code\\cbem_new\\3d.dat") + ".plt");
		return 1;
	}
}
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( MathTest, "alltests" );
//CPPUNIT_TEST_SUITE_REGISTRATION( MeshTest );


