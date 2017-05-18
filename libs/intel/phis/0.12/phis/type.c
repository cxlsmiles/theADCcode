#include <stdio.h>
#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _
#define phis_init F77_FUNC_(phis_init,PHIS_INIT)

using namespace std;

extern "C" {
struct SParamBase {
 int iVersion;
};
struct SParamSet1: public SParamBase {
 int flags;
};
struct SParamSet2: SParamSet1 {
 char *path;
};
}

extern"C" {
void phis_init(SParamBase *par){
 switch(par->iVersion){
  case 0:{
    printf("1 para\n");
      SParamSet1 *set1 = 0;
      set1 = (SParamSet1*)par;
      printf("flags[0]=%d\n", set1->flags);
    }
    break;
  case 1: {
      SParamSet2 *set2 = 0;
      printf("2 para\n");
      set2 = (SParamSet2*)par;
      printf("flags[0]=%d, par2=%s\n", set2->flags, set2->path);
    }
    break;
  default:
    printf("error: version not supported\n");
    break;
 }
}
}

/*int main(){
 int flags1[] = {1,2,3};
 int flags2[] = {6,5,4};
 char path[] = "C:\\path";

 SParamSet1 par1;
 par1.iVersion = 0;
 par1.flags = flags1;
 phis_init(&par1);

 SParamSet2 par2;
 par2.iVersion = 1;
 par2.flags = flags2;
 par2.path = path;
 phis_init(&par2);}
 */

