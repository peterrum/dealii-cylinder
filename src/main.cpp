#include <deal.II/base/point.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/tria.h>
#include <deal.II/numerics/data_out.h>

#include <stdio.h>
#include <vector>

using namespace dealii;

const double r1 = 1.0;
const double r2 = 0.8;
const double l  = 5.0;

void
create_reference_cylinder(const bool                 do_transition,
                          const unsigned int         n_sections,
                          std::vector<Point<3>> &    vertices_3d,
                          std::vector<CellData<3>> & cell_data_3d)
{
  if(do_transition)
    printf("WARNING: Transition has not been implemented yet (TODO)!\n");

  // position of auxiliary point to achieve an angle of 120 degrees in corner
  // of inner cell
  const double radius = 1;
  const double ycord  = 0.55 * radius * std::cos(numbers::PI / 12) /
                       (std::sin(numbers::PI / 12) + std::cos(numbers::PI / 12));
  // vertices for quarter of circle
  std::vector<Point<2>> vertices{{0, 0},
                                 {0.55 * radius, 0},
                                 {ycord, ycord},
                                 {radius, 0},
                                 {radius * std::sqrt(0.5), radius * std::sqrt(0.5)}};

  // create additional vertices for other three quarters of circle -> gives 17
  // vertices in total
  for(unsigned int a = 1; a < 4; ++a)
  {
    Tensor<2, 2> transform;
    transform[0][0] = a == 2 ? -1. : 0;
    transform[1][0] = a == 2 ? 0 : (a == 1 ? 1 : -1);
    transform[0][1] = -transform[1][0];
    transform[1][1] = transform[0][0];
    for(unsigned int i = 1; i < 5; ++i)
      vertices.push_back(Point<2>(transform * vertices[i]));
  }

  // create 12 cells for 2d mesh on base; the first four elements are at the
  // center of the circle
  std::vector<CellData<2>> cell_data(12);
  cell_data[0].vertices[0] = 0;
  cell_data[0].vertices[1] = 1;
  cell_data[0].vertices[2] = 5;
  cell_data[0].vertices[3] = 2;
  cell_data[1].vertices[0] = 9;
  cell_data[1].vertices[1] = 0;
  cell_data[1].vertices[2] = 6;
  cell_data[1].vertices[3] = 5;
  cell_data[2].vertices[0] = 10;
  cell_data[2].vertices[1] = 13;
  cell_data[2].vertices[2] = 9;
  cell_data[2].vertices[3] = 0;
  cell_data[3].vertices[0] = 13;
  cell_data[3].vertices[1] = 14;
  cell_data[3].vertices[2] = 0;
  cell_data[3].vertices[3] = 1;

  // the next 8 elements describe the rim; we take one quarter of the circle
  // in each loop iteration
  for(unsigned int a = 0; a < 4; ++a)
  {
    cell_data[4 + a * 2].vertices[0] = 1 + a * 4;
    cell_data[4 + a * 2].vertices[1] = 3 + a * 4;
    cell_data[4 + a * 2].vertices[2] = 2 + a * 4;
    cell_data[4 + a * 2].vertices[3] = 4 + a * 4;
    cell_data[5 + a * 2].vertices[0] = 2 + a * 4;
    cell_data[5 + a * 2].vertices[1] = 4 + a * 4;
    AssertIndexRange(4 + a * 4, vertices.size());
    cell_data[5 + a * 2].vertices[2] = a == 3 ? 1 : 5 + a * 4;
    cell_data[5 + a * 2].vertices[3] = a == 3 ? 3 : 7 + a * 4;
  }
  SubCellData subcell_data;
  GridReordering<2>::reorder_cells(cell_data, true);

  Triangulation<2> tria_2d;
  tria_2d.create_triangulation(vertices, cell_data, subcell_data);

  vertices_3d.clear();
  vertices_3d.resize((n_sections + 1) * tria_2d.n_vertices());
  cell_data_3d.clear();
  cell_data_3d.resize(n_sections * cell_data.size());

  for(unsigned int s = 0; s <= n_sections; s++)
  {
    const double       beta  = (1.0 * s) / n_sections;
    const unsigned int shift = s * tria_2d.n_vertices();
    for(unsigned int i = 0; i < tria_2d.n_vertices(); ++i)
    {
      vertices_3d[shift + i][0] = tria_2d.get_vertices()[i][0];
      vertices_3d[shift + i][1] = tria_2d.get_vertices()[i][1];
      vertices_3d[shift + i][2] = beta;
    }
  }


  for(unsigned int s = 0; s < n_sections; s++)
    for(unsigned int i = 0; i < cell_data.size(); ++i)
    {
      for(unsigned int v = 0; v < 4; ++v)
        cell_data_3d[cell_data.size() * s + i].vertices[v + 0] =
          (s + 0) * vertices.size() + cell_data[i].vertices[v];
      for(unsigned int v = 0; v < 4; ++v)
        cell_data_3d[cell_data.size() * s + i].vertices[4 + v] =
          (s + 1) * vertices.size() + cell_data[i].vertices[v];
    }
}

void
create_cylinder(const double               radius1,
                const double               radius2,
                const double               length,
                const bool                 do_transition,
                std::vector<Point<3>> &    vertices_3d,
                std::vector<CellData<3>> & cell_data_3d)
{
  // create reference cylinder with n_sections subdivisions
  int n_sections = length / std::min(radius1, radius2);
  create_reference_cylinder(do_transition, n_sections, vertices_3d, cell_data_3d);

  // transform cylinder (here: simple scaling)
  for(auto & point : vertices_3d)
  {
    const double beta  = point[2];
    const double alpha = 1.0 - beta;

    const double r = beta * radius2 + alpha * radius1;
    const double l = beta * length;

    point[0] = r * point[0];
    point[1] = r * point[1];
    point[2] = l;
  }
}

int
main(int argc, char ** argv)
{
  // default size of cylinder
  double radius_1 = r1, radius_2 = r2, length = l;
  bool   do_transition = false;

  // read new sizes of geometry from command line
  if(argc > 1)
    radius_1 = atof(argv[1]);
  if(argc > 2)
    radius_2 = atof(argv[2]);
  if(argc > 3)
    length = atof(argv[3]);
  if(argc > 4)
    do_transition = atoi(argv[4]);

  std::vector<CellData<3>> cell_data_3d;
  std::vector<Point<3>>    vertices_3d;
  create_cylinder(radius_1, radius_2, length, do_transition, vertices_3d, cell_data_3d);

  // now actually create the triangulation
  SubCellData      subcell_data;
  Triangulation<3> tria;

  try
  {
    tria.create_triangulation(vertices_3d, cell_data_3d, subcell_data);
  }
  catch(const std::exception & e)
  {
    AssertThrow(false, ExcMessage(e.what()));
  }

  DataOut<3> data_out_3d;
  data_out_3d.attach_triangulation(tria);
  data_out_3d.build_patches();
  std::ofstream file_3d("grid_3d.vtk");
  data_out_3d.write_vtk(file_3d);
  tria.refine_global();

  // output triangulation to VTK-file (Paraview)
}
