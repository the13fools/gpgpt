
#include <iostream>

template<int M>
class Impl
{

    public: 


    template<int K>
    void functest2()
    {
        std::cout << "K " << K << " M " << M << std::endl;
    } 

};

template<int N> 
class Top
{
    public:
    
    void functest(Impl<N>& tt )
    {

        // OMG THE WORST!!!!
        tt.template functest2<3>();    
    }
};

int main()
{
    Top<5> t2;
    Impl<5> tt;
    t2.functest(tt);

    return 0;
}