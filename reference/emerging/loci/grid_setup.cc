#include <iostream>
using namespace std ;

main() {
  int edgeNodes = 46 ;
  cout << "1" << endl ;
  cout << edgeNodes << ' ' << edgeNodes << ' ' << edgeNodes << endl ;
  cout.precision(16) ;
  for(int i=0;i<edgeNodes;++i)
    for(int j=0;j<edgeNodes;++j)
      for(int k=0;k<edgeNodes;++k) {
        cout << 1.125*double(k)/double(edgeNodes-1) << endl ;
      }
  for(int i=0;i<edgeNodes;++i)
    for(int j=0;j<edgeNodes;++j)
      for(int k=0;k<edgeNodes;++k) {
        cout << 1.125*double(j)/double(edgeNodes-1) << endl ;
      }
  for(int i=0;i<edgeNodes;++i)
    for(int j=0;j<edgeNodes;++j)
      for(int k=0;k<edgeNodes;++k) {
        cout << 1.125*double(i)/double(edgeNodes-1) << endl ;
      }
}
