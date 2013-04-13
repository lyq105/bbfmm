// Using cppunit test the mesh io class.
// co
//
#include "cppunit/extensions/HelperMacros.h"

class MeshTest: public CppUnit::TestFixture
{

	CPPUNIT_TEST_SUITE( MeshTest );
	
	
	CPPUNIT_TEST( test_read );
	
	
	CPPUNIT_TEST_SUITE_END();


	public:
	MeshTest(){}
	// 初始化类函数
	void setUp();
	// 后处理函数
	void tearDown();
  
	// 测试用例
	void test_read();
};

CPPUNIT_TEST_SUITE_NAMED_REGISTRATION( MathTest, "alltests" );
//CPPUNIT_TEST_SUITE_REGISTRATION( MeshTest );


