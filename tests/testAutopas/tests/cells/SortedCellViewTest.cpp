/**
 * @file SortedCellViewTest.cpp
 * @author C. Menges
 * @date 26.05.2019
 */

#include "SortedCellViewTest.h"

#include "autopas/cells/FullParticleCell.h"
#include "autopas/cells/SortedCellView.h"
#include "testingHelpers/commonTypedefs.h"

TEST_F(SortedCellViewTest, testParticleAccess) {
  auto fpc = autopas::FullParticleCell<Particle>();
  Particle p1 = Particle();
  fpc.addParticle(p1);
  auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {1., 0., 0.});
  EXPECT_EQ(fspc._particles.size(), 1);
  std::array<double, 3> force{3.1416, 2.7183, 9.8067};
  fspc._particles.front().second->addF(force);
  EXPECT_THAT(fpc._particles.front().getF(), testing::ContainerEq(force));
}

TEST_F(SortedCellViewTest, testParticleSorting) {
  auto fpc = autopas::FullParticleCell<Particle>();
  Particle p1 = Particle({0., 3., 0.}, {0., 0., 0.}, 0);
  fpc.addParticle(p1);
  Particle p2 = Particle({2., 1., 2.}, {0., 0., 0.}, 2);
  fpc.addParticle(p2);
  Particle p3 = Particle({1., 2., 1.}, {0., 0., 0.}, 1);
  fpc.addParticle(p3);
  Particle p4 = Particle({3., 0., 3.}, {0., 0., 0.}, 3);
  fpc.addParticle(p4);

  unsigned int id = 0u;

  {
    auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {1., 0., 0.});
    EXPECT_EQ(fspc.numParticles(), 4);

    for (auto &p : fspc._particles) {
      EXPECT_DOUBLE_EQ(p.first, static_cast<double>(id));
      EXPECT_EQ(p.second->getID(), id);
      id++;
    }
  }
  {
    auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, {0., 1., 0.});
    EXPECT_EQ(fspc.numParticles(), 4);

    for (auto &p : fspc._particles) {
      id--;
      EXPECT_DOUBLE_EQ(p.first, static_cast<double>(3u - id));
      EXPECT_EQ(p.second->getID(), id);
    }
  }
}
static double rand_doub(double fmin, double fmax){
  double f = (double) rand()/((double)RAND_MAX);
  return fmin + f*(fmax - fmin);
}
static std::string arrtostr(std::array<double, 3> a){
  return "{"+std::to_string(a[0])+", "+std::to_string(a[1])+", "+std::to_string(a[2])+"}";
}
/**
 * Tests if iterating only until projected values from sorting have more distance than cutoff
 * will not ignore particle pairs that have distance leq cutoff.
 * Tested with 10,000 particles, but reduced the max_parts to 1000 so that this doesn't take as much time.
 */
TEST_F(SortedCellViewTest, testCutoff){
  const int max_parts = 1000;
  for(int num_parts : {3,20, 125, 300, 1000}){
    for(double cutoff : {0.5, 1.0, 2.0}){
      auto fpc = autopas::FullParticleCell<Particle>();
      std::array<Particle, max_parts> particles;
      for(int i = 0; i < num_parts; i++){
        double xr = rand_doub(0.0, 4.0);
        double yr = rand_doub(0.0, 4.0);
        double zr = rand_doub(0.0, 4.0);
        particles[i] = Particle({xr,yr,zr},{0.,0.,0.},i);
        fpc.addParticle(particles[i]);
      }
      {
        //auto cutoff = 1.0;
        auto fspc = autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>>(fpc, autopas::utils::ArrayMath::normalize(std::array<double,3>{1., 1., 1.}));
        //std::array<std::array<bool, max_parts>, max_parts> visited{false};
        //std::array<std::array<bool, max_parts>, max_parts> shouldVisit{false};
        //std::array<std::array<bool, max_parts>, max_parts> sq_visit{false};
        auto outer = fspc._particles.begin();
        for (; outer != fspc._particles.end(); ++outer) {
          Particle &p1 = *outer->second;

          auto inner = outer;
          ++inner;
          bool skip = false;
          for (; inner != fspc._particles.end(); ++inner) {
            Particle &p2 = *inner->second;
            if (std::abs(outer->first - inner->first) > cutoff) {
              skip = true;
            }
            //visited[outer->second->getID()][inner->second->getID()] = !skip;
            std::array<double, 3> dist = autopas::utils::ArrayMath::sub(inner->second->getR(), outer->second->getR());
            double squared_norm = autopas::utils::ArrayMath::dot(dist, dist);
            double _norm = autopas::utils::ArrayMath::L2Norm(dist);
            //shouldVisit[outer->second->getID()][inner->second->getID()] = _norm <= cutoff;
            //sq_visit[outer->second->getID()][inner->second->getID()] = squared_norm <= cutoff;
            if(skip){
              //EXPECT_FALSE(shouldVisit[outer->second->getID()][inner->second->getID()]);
              //EXPECT_FALSE(sq_visit[outer->second->getID()][inner->second->getID()]);
              EXPECT_FALSE(_norm <= cutoff);
              //For Debugging:
              /*EXPECT_TRUE(false) << "Cutoff " << cutoff << std::endl
                              << "Proj_D " << std::abs(outer->first - inner->first) << std::endl
                              << "Norm   " << _norm << std::endl
                              << "PosI   " << arrtostr(inner->second->getR()) << std::endl
                              << "PosO   " << arrtostr(outer->second->getR());*/
            }
            EXPECT_LE(std::abs(outer->first - inner->first), _norm) << std::endl << "Cutoff " << cutoff;
          }
        }
      }
    }
  }

}
static int get_from_3d(std::array<int, 3> d3, std::array<int, 3> dims){
  return d3[2] + dims[2]*(d3[1] + dims[1]*(d3[0]));
}
static void checkCutoffs(autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>> a, autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>> b){
  for(double cutoff : {/*0.5, */1.0, 2.0}) {
    //Copied from CellFunctor (to really check real behaviour):
    for (auto &outer : a._particles) {
      Particle &p1 = *outer.second;
      bool skip = false;
      for (auto &inner : b._particles) {
        if (std::abs(outer.first - inner.first) > cutoff) {
          skip = true;
        }
        Particle &p2 = *inner.second;
        std::array<double, 3> dist = autopas::utils::ArrayMath::sub(p2.getR(), p1.getR());
        double _norm = autopas::utils::ArrayMath::L2Norm(dist);
        if(skip){
          EXPECT_GT(_norm, cutoff);
        }
        EXPECT_LE(std::abs(outer.first - inner.first), _norm);
      }
    }
  }
}
TEST_F(SortedCellViewTest, testCutoffCellPair){
  const int max_parts = 1000;
  int num_parts = max_parts;
  std::array<int, 3> cell_dims = {5, 5, 5};
  const int cell_amount = 5*5*5;
  EXPECT_EQ(cell_amount, cell_dims[0]*cell_dims[1]*cell_dims[2]);
  std::array<autopas::FullParticleCell<Particle>,cell_amount> fpcs{autopas::FullParticleCell<Particle>()};
  std::array<Particle, max_parts> particles;
  for(int i = 0; i < num_parts; i++){
    double xr = rand_doub(0.0, (double)cell_dims[0]);
    double yr = rand_doub(0.0, (double)cell_dims[1]);
    double zr = rand_doub(0.0, (double)cell_dims[2]);
    particles[i] = Particle({xr,yr,zr},{0.,0.,0.},i);
    fpcs[get_from_3d({(int)xr, (int)yr, (int)zr}, cell_dims)].addParticle(particles[i]);
  }
  double cutoff = 1.;
  for(int xi = 0, ti = 0; xi < cell_dims[0]; ++xi){
    for (int yi = 0; yi < cell_dims[1]; ++yi) {
      for (int zi = 0; zi < cell_dims[2]; ++zi, ++ti) {
        EXPECT_LT(ti, cell_amount);
        for(int xi2 = xi, ti2 = ti+1; xi2 < cell_dims[0]; ++xi2){
          for(int yi2 = yi; yi2 < cell_dims[1]; ++yi2){
            for(int zi2 = zi + 1; zi2 < cell_dims[2]; ++zi2){
              std::array<double, 3> dist = {(double) xi2-xi, (double) yi2-yi, (double) zi2-zi};
              std::array<double, 3> normalized = autopas::utils::ArrayMath::normalize(dist);
              autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>> scv1(fpcs[ti], normalized);
              //std::array<double, 3> neg_normal = autopas::utils::ArrayMath::mulScalar(normalized, -1.);
              autopas::SortedCellView<Particle, autopas::FullParticleCell<Particle>> scv2(fpcs[ti2], normalized);
              checkCutoffs(scv1, scv2);
            }
          }
        }
      }
    }
  }
}