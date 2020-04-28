#include <math.h>
#include <vector>
#include <fstream>

#include <cxxopts.hpp>
#include <spdlog/spdlog.h>

#include <xtensor/xio.hpp>

#include <mure/mure.hpp>
#include "coarsening.hpp"
#include "refinement.hpp"
#include "criteria.hpp"

#include <chrono>


/// Timer used in tic & toc
auto tic_timer = std::chrono::high_resolution_clock::now();

/// Launching the timer
void tic()
{
    tic_timer = std::chrono::high_resolution_clock::now();
}


/// Stopping the timer and returning the duration in seconds
double toc()
{
    const auto toc_timer = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_span = toc_timer - tic_timer;
    return time_span.count();
}

double exact_solution(double x, double t)   {

    double u = 0;

    // { // Hyperbolic tangent

    //     double sigma = 20.0;

    //     if (t <= 0.0)
    //         return 0.5 * (1.0 + tanh(sigma * x));
    //     else
    //     {   // We proceed by dicothomy
    //         double a = -3.2;
    //         double b =  3.2;

    //         double tol = 1.0e-8;

    //         auto F = [sigma, x, t] (double y)   {
    //             return y + 0.5 * (1.0 + tanh(sigma * y))*t - x;
    //         };
    //         double res = 0.0;

    //         while (b-a > tol)   {
    //             double mean = 0.5 * (b + a);
    //             double eval = F(mean);
    //             if (eval <= 0.0)
    //                 a = mean;
    //             else
    //                 b = mean;
    //             res = mean;
    //         }

    //         return 0.5 * (1.0 + tanh(sigma * res));

    //     }
        
    // }

    // double rhoL = 1.0;
    // double rhoR = 0.0;
    // double x0 = 0.0;

    //double vshock = 0.5 * (rhoL + rhoR);

    //return ((x-0.5*t) <= x0) ? rhoL : rhoR;
    //return ((x-0.5*t) <= -1.0) ? 0 : ((x-0.5*t) <= 1.0 ? 1.0 : 0.0);

    //return (x <= x0 + vshock * t) ? rhoL : rhoR;
    // {
    //     double x0L = -0.2;
    //     double x0R =  0.2;
    //     return ((x - 0.5*t) < x0L) ? 0.0 : (((x - 0.5*t) < x0R ? 1.0 : 0.0));
    // }

    // double sigma = 0.5;
    // double rhoL = 0.0;
    // double rhoC = 1.0;
    // double rhoR = 0.0;

    // return (x + sigma <= rhoL * t) ? rhoL : ((x + sigma <= rhoC*t) ? (x+sigma)/t : ((x-sigma <= t/2*(rhoC + rhoR)) ? rhoC : rhoR ));


    // // x = x - 0.5 * t;
    // // t = 0.0;

    if (x >= -1 and x < t)
    {
        u = (1 + x) / (1 + t);
    }
    
    if (x >= t and x < 1)
    {
        u = (1 - x) / (1 - t);
    }

    //u = exp(-20.0 * (x-0.75*t) * (x-0.75*t));

    // u = 0.0;
    
    return u;
}

double flux(double u)   {
    return 0.5 * u * u;
    //return 0.75 * u;
}

template<class Config>
auto init_f(mure::Mesh<Config> &mesh, double t)
{
    constexpr std::size_t nvel = 2;
    mure::BC<1> bc{ {{ {mure::BCType::neumann, 0.0},
                       {mure::BCType::neumann, 0.0},
                    }} };

    mure::Field<Config, double, nvel> f("f", mesh, bc);
    f.array().fill(0);

    mesh.for_each_cell([&](auto &cell) {
        auto center = cell.center();
        auto x = center[0];
        double u = 0;

        // if (x >= -1 and x < t)
        // {
        //     u = (1 + x) / (1 + t);
        // }
        // if (x >= t and x < 1)
        // {
        //     u = (1 - x) / (1 - t);
        // }

        u = exact_solution(x, 0.0);

        //double u = exp(-20.0 * x * x);

        double v = flux(u);//.5 * u; 
        //double v = .5 * u * u;

        f[cell][0] = .5 * (u + v);
        f[cell][1] = .5 * (u - v);
    });

    return f;
}

template<class Field, class interval_t, class FieldTag>
xt::xtensor<double, 1> prediction(const Field& f, std::size_t level_g, std::size_t level, const interval_t &i, const std::size_t item, 
                                  const FieldTag & tag, std::map<std::tuple<std::size_t, std::size_t, std::size_t, interval_t>, 
                                  xt::xtensor<double, 1>> & mem_map)
{

    // We check if the element is already in the map
    auto it = mem_map.find({item, level_g, level, i});
    if (it != mem_map.end())   {
        //std::cout<<std::endl<<"Found by memoization";
        return it->second;
    }
    else {

        auto mesh = f.mesh();
        xt::xtensor<double, 1> out = xt::empty<double>({i.size()/i.step});//xt::eval(f(item, level_g, i));
        auto mask = mesh.exists(level_g + level, i);

        // std::cout << level_g + level << " " << i << " " << mask << "\n"; 
        if (xt::all(mask))
        {         
            return xt::eval(f(item, level_g + level, i));
        }

        auto step = i.step;
        auto ig = i / 2;
        ig.step = step >> 1;
        xt::xtensor<double, 1> d = xt::empty<double>({i.size()/i.step});

        for (int ii=i.start, iii=0; ii<i.end; ii+=i.step, ++iii)
        {
            d[iii] = (ii & 1)? -1.: 1.;
        }

    
        auto val = xt::eval(prediction(f, level_g, level-1, ig, item, tag, mem_map) - 1./8 * d * (prediction(f, level_g, level-1, ig+1, item, tag, mem_map) 
                                                                                       - prediction(f, level_g, level-1, ig-1, item, tag, mem_map)));
        

        xt::masked_view(out, !mask) = xt::masked_view(val, !mask);
        for(int i_mask=0, i_int=i.start; i_int<i.end; ++i_mask, i_int+=i.step)
        {
            if (mask[i_mask])
            {
                out[i_mask] = f(item, level_g + level, {i_int, i_int + 1})[0];
            }
        }

        // The value should be added to the memoization map before returning
        return mem_map[{item, level_g, level, i}] = out;

        //return out;
    }

}


template<class Field, class interval_t>
xt::xtensor<double, 2> prediction_all(const Field& f, std::size_t level_g, std::size_t level, const interval_t &i, 
                                  std::map<std::tuple<std::size_t, std::size_t, interval_t>, 
                                  xt::xtensor<double, 2>> & mem_map)
{

    using namespace xt::placeholders;
    // We check if the element is already in the map
    auto it = mem_map.find({level_g, level, i});
    if (it != mem_map.end())
    {
        return it->second;
    }
    else
    {
        auto mesh = f.mesh();
        std::vector<std::size_t> shape = {i.size(), 2};
        xt::xtensor<double, 2> out = xt::empty<double>(shape);
        auto mask = mesh.exists(level_g + level, i);

        xt::xtensor<double, 2> mask_all = xt::empty<double>(shape);
        xt::view(mask_all, xt::all(), 0) = mask;
        xt::view(mask_all, xt::all(), 1) = mask;

        if (xt::all(mask))
        {         
            return xt::eval(f(level_g + level, i));
        }

        auto ig = i >> 1;
        ig.step = 1;

        xt::xtensor<double, 2> val = xt::empty<double>(shape);
        auto current = xt::eval(prediction_all(f, level_g, level-1, ig, mem_map));
        auto left = xt::eval(prediction_all(f, level_g, level-1, ig-1, mem_map));
        auto right = xt::eval(prediction_all(f, level_g, level-1, ig+1, mem_map));

        std::size_t start_even = (i.start&1)? 1: 0;
        std::size_t start_odd = (i.start&1)? 0: 1;
        std::size_t end_even = (i.end&1)? ig.size(): ig.size()-1;
        std::size_t end_odd = (i.end&1)? ig.size()-1: ig.size();
        xt::view(val, xt::range(start_even, _, 2)) = xt::view(current - 1./8 * (right - left), xt::range(start_even, _));
        xt::view(val, xt::range(start_odd, _, 2)) = xt::view(current + 1./8 * (right - left), xt::range(_, end_odd));

        xt::masked_view(out, !mask_all) = xt::masked_view(val, !mask_all);
        for(int i_mask=0, i_int=i.start; i_int<i.end; ++i_mask, ++i_int)
        {
            if (mask[i_mask])
            {
                xt::view(out, i_mask) = xt::view(f(level_g + level, {i_int, i_int + 1}), 0);
            }
        }

        // The value should be added to the memoization map before returning
        return out;// mem_map[{level_g, level, i, ig}] = out;
    }
}

template<class Field, class FieldTag>
void one_time_step(Field &f, const FieldTag & tag, double s)
{
    constexpr std::size_t nvel = Field::size;
    double lambda = 1.;//, s = 1.0;
    auto mesh = f.mesh();
    auto max_level = mesh.max_level();

    mure::mr_projection(f);
    f.update_bc();
    mure::mr_prediction(f);


    // MEMOIZATION
    // All is ready to do a little bit  of mem...
    using interval_t = typename Field::Config::interval_t;
    std::map<std::tuple<std::size_t, std::size_t, std::size_t, interval_t>, xt::xtensor<double, 1>> memoization_map;
    memoization_map.clear(); // Just to be sure...

    Field new_f{"new_f", mesh};
    new_f.array().fill(0.);

    for (std::size_t level = 0; level <= max_level; ++level)
    {
        auto exp = mure::intersection(mesh[mure::MeshType::cells][level],
                                      mesh[mure::MeshType::cells][level]);
        exp([&](auto, auto &interval, auto) {
            auto i = interval[0];


            // STREAM

            std::size_t j = max_level - level;

            double coeff = 1. / (1 << j);

            // This is the STANDARD FLUX EVALUATION
            
            auto fp = f(0, level, i) + coeff * (prediction(f, level, j, i*(1<<j)-1, 0, tag, memoization_map)
                                             -  prediction(f, level, j, (i+1)*(1<<j)-1, 0, tag, memoization_map));

            auto fm = f(1, level, i) - coeff * (prediction(f, level, j, i*(1<<j), 1, tag, memoization_map)
                                             -  prediction(f, level, j, (i+1)*(1<<j), 1, tag, memoization_map));
            
            
            // This is the CHEAP FLUX EVALUATION
            
            // auto fp = f(0, level, i); // Just to give the shape ....
            // auto fm = f(1, level, i);

            // auto exist_0_m1 = mesh.exists(level, i - 1, mure::MeshType::cells);
            // auto exist_0_p1 = mesh.exists(level, i + 1, mure::MeshType::cells);

            // auto exist_up_m1 = mesh.exists(level + 1, 2*i - 1, mure::MeshType::cells);
            // auto exist_up_p1 = mesh.exists(level + 1, 2*i + 1, mure::MeshType::cells);

            // auto exist_down_m1 = mesh.exists(level - 1, i/2 - 1, mure::MeshType::cells); // Verify this division by 2... sometimes is problematic
            // auto exist_down_p1 = mesh.exists(level - 1, i/2 + 1, mure::MeshType::cells);

            // // The left neigh at the same level exists
            // xt::masked_view(fp, exist_0_m1) = (1.0 - coeff) * f(0, level, i) + coeff * f(0, level, i - 1);
            // // The right neigh at the same level exists
            // xt::masked_view(fm, exist_0_p1) = (1.0 - coeff) * f(1, level, i) + coeff * f(1, level, i + 1);


            // // This is problematic... ASK
            // xt::masked_view(xt::masked_view(fp, !exist_0_m1), 
            //                                     exist_up_m1) = (1.0 - coeff) * f(0, level, i) + coeff * f(0, level + 1, 2*i - 1); 
                                                
            // xt::masked_view(xt::masked_view(fm, !exist_0_p1), 
            //                                         exist_up_p1) = (1.0 - coeff) * f(1, level, i);// + coeff * f(1, level + 1, 2*i + 1); 


            


            // COLLISION    

            auto uu = xt::eval(fp + fm);
            auto vv = xt::eval(lambda * (fp - fm));

            
            //vv = (1 - s) * vv + s * 0.75 * uu;

            vv = (1 - s) * vv + s * .5 * uu * uu;

            new_f(0, level, i) = .5 * (uu + 1. / lambda * vv);
            new_f(1, level, i) = .5 * (uu - 1. / lambda * vv);
        });
    }

    std::swap(f.array(), new_f.array());
}

template<class Field>
void save_solution(Field &f, double eps, std::size_t ite, std::string ext)
{
    using Config = typename Field::Config;
    auto mesh = f.mesh();
    std::size_t min_level = mesh.min_level();
    std::size_t max_level = mesh.max_level();

    std::stringstream str;
    str << "LBM_D1Q2_Burgers_" << ext << "_lmin_" << min_level << "_lmax-" << max_level << "_eps-"
        << eps << "_ite-" << ite;

    auto h5file = mure::Hdf5(str.str().data());
    h5file.add_mesh(mesh);
    mure::Field<Config> level_{"level", mesh};
    mure::Field<Config> u{"u", mesh};
    mesh.for_each_cell([&](auto &cell) {
        level_[cell] = static_cast<double>(cell.level);
        u[cell] = f[cell][0] + f[cell][1];
    });
    h5file.add_field(u);
    h5file.add_field(f);
    h5file.add_field(level_);
}

template<class Field>
void save_refined_solution(Field &f, std::size_t min_level, std::size_t max_level, double eps, std::size_t ite, std::string ext="")
{
    using Config = typename Field::Config;
    auto mesh = f.mesh();

    std::stringstream str;
    str << "LBM_D1Q2_Burgers_refined_solution_" << ext << "_lmin_" << min_level << "_lmax-" << max_level << "_eps-"
        << eps << "_ite-" << ite;

    auto h5file = mure::Hdf5(str.str().data());
    h5file.add_mesh(mesh);
    mure::Field<Config> level_{"level", mesh};
    mure::Field<Config> u{"u", mesh};
    mesh.for_each_cell([&](auto &cell) {
        level_[cell] = static_cast<double>(cell.level);
        u[cell] = f[cell][0] + f[cell][1];
    });
    h5file.add_field(u);
    h5file.add_field(f);
    h5file.add_field(level_);
}

// template<class Field, class FieldTag>
// double compute_error(const Field & f, const FieldTag & tag, double t)
// {
//     double error_to_return = 0.0;

//     // Getting ready for memoization
//     using interval_t = typename Field::Config::interval_t;
//     std::map<std::tuple<std::size_t, std::size_t, std::size_t, interval_t>, xt::xtensor<double, 1>> memoization_map;
//     memoization_map.clear();


//     auto mesh = f.mesh();
//     auto max_level = mesh.max_level();

//     double dx = 1 << max_level;


//     Field error_cell_by_cell{"error_cell_by_cell", mesh};
//     error_cell_by_cell.array().fill(0.);


//     auto subset = intersection(mesh.initial_mesh(), mesh.initial_mesh()).on(max_level);
//     subset([&](auto, auto &interval, auto) {
//         auto i = interval[0];


//         std::cout<<"\n\nHere  "<<i;

//         auto fp = prediction(f, max_level, 0, i, 0, tag, memoization_map); // BUG ICI
//         //auto fm = prediction(f, max_level, 0, i, 1, tag, memoization_map);

//         // CELA IL FAUT LE TRAITER MAIS ON VERRA APRES...

//         // auto rho = xt::eval(fp + fm);

//         // error_cell_by_cell(0, max_level, i) = dx * xt::abs(rho);
//     });

//     return error_to_return;
    
    
//     //return xt::sum(xt::view(error_cell_by_cell, 0, max_level, xt::all));
// }

// template<class Field, class FieldR>
template<class Config, class FieldR>
std::array<double, 2> compute_error(mure::Field<Config, double, 2> &f, FieldR & fR, double t)
{

    auto mesh = f.mesh();

    auto meshR = fR.mesh();
    auto max_level = meshR.max_level();

    fR.update_bc();    

    mure::mr_projection(f);    
    f.update_bc(); // Important especially when we enforce Neumann...for the Riemann problem
    mure::mr_prediction(f);  // C'est supercrucial de le faire.


    // Getting ready for memoization
    // using interval_t = typename Field::Config::interval_t;
    using interval_t = typename Config::interval_t;
    std::map<std::tuple<std::size_t, std::size_t, interval_t>, xt::xtensor<double, 2>> error_memoization_map;
    error_memoization_map.clear();

    double error = 0; // To return
    double diff = 0.0;


    double dx = 1.0 / (1 << max_level);

    for (std::size_t level = 0; level <= max_level; ++level)
    {
        auto exp = mure::intersection(meshR[mure::MeshType::cells][max_level],
                                      mesh[mure::MeshType::cells][level])
                  .on(max_level);

        exp([&](auto, auto &interval, auto) {
            auto i = interval[0];
            auto j = max_level - level;

            auto sol  = prediction_all(f, level, j, i, error_memoization_map);
            auto solR = xt::view(fR(max_level, i), xt::all(), xt::range(0, 2));


            xt::xtensor<double, 1> x = dx*xt::linspace<int>(i.start, i.end - 1, i.size()) + 0.5*dx;
            // xt::xtensor<double, 1> uexact = (x >= -1.0 and x < t) * ((1 + x) / (1 + t)) + 
            //                                 (x >= t and x < 1) * (1 - x) / (1 - t);

            xt::xtensor<double, 1> uexact = xt::zeros<double>(x.shape());

            for (std::size_t idx = 0; idx < x.shape()[0]; ++idx)    {
                uexact[idx] = exact_solution(x[idx], t); // We can probably do better
            }

            error += xt::sum(xt::abs(xt::flatten(xt::view(fR(max_level, i), xt::all(), xt::range(0, 1)) + xt::view(fR(max_level, i), xt::all(), xt::range(1, 2))) 
                             - uexact))[0];


            diff += xt::sum(xt::abs(xt::flatten(xt::view(sol, xt::all(), xt::range(0, 1)) + xt::view(sol, xt::all(), xt::range(1, 2))) - xt::flatten(xt::view(fR(max_level, i), xt::all(), xt::range(0, 1)) + xt::view(fR(max_level, i), xt::all(), xt::range(1, 2)))))[0];


        });
    }

    return {dx * error, dx * diff}; // Normalization by dx before returning
    // I think it is better to do the normalization at the very end ... especially for round-offs    
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("lbm_d1q2_burgers",
                             "Multi resolution for a D1Q2 LBM scheme for Burgers equation");

    options.add_options()
                       ("min_level", "minimum level", cxxopts::value<std::size_t>()->default_value("2"))
                       ("max_level", "maximum level", cxxopts::value<std::size_t>()->default_value("10"))
                       ("epsilon", "maximum level", cxxopts::value<double>()->default_value("0.01"))
                       ("s", "relaxation parameter", cxxopts::value<double>()->default_value("1.0"))
                       ("log", "log level", cxxopts::value<std::string>()->default_value("warning"))
                       ("h, help", "Help");

    try
    {
        auto result = options.parse(argc, argv);

        if (result.count("help"))
            std::cout << options.help() << "\n";
        else
        {
            std::map<std::string, spdlog::level::level_enum> log_level{{"debug", spdlog::level::debug},
                                                               {"warning", spdlog::level::warn}};
            constexpr size_t dim = 1;
            using Config = mure::MRConfig<dim, 2>;

            spdlog::set_level(log_level[result["log"].as<std::string>()]);
            std::size_t min_level = result["min_level"].as<std::size_t>();
            std::size_t max_level = result["max_level"].as<std::size_t>();
            double eps = result["epsilon"].as<double>();
            double s = result["s"].as<double>();


            mure::Box<double, dim> box({-3}, {3});
            mure::Mesh<Config> mesh{box, min_level, max_level};
            mure::Mesh<Config> meshR{box, max_level, max_level}; // This is the reference scheme

            // for (std::size_t level = min_level; level <= max_level; ++level) {

            //     auto cells = intersection(mesh[mure::MeshType::cells][level], 
            //                               mesh[mure::MeshType::cells][level]);

            //     auto allcells = intersection(mesh[mure::MeshType::all_cells][level], 
            //                                  mesh[mure::MeshType::all_cells][level]);
            //     auto projcells = intersection(mesh[mure::MeshType::proj_cells][level], 
            //                                   mesh[mure::MeshType::proj_cells][level]);

            //     auto cellsandghosts = intersection(mesh[mure::MeshType::cells_and_ghosts][level], 
            //                                        mesh[mure::MeshType::cells_and_ghosts][level]);


            //     cells([&](auto, auto &interval, auto) {
            //         auto i = interval[0];
            //         std::cout<<std::endl<<"Level "<<level<<"Cells "<<i;
            //     });
                
            //     allcells([&](auto, auto &interval, auto) {
            //         auto i = interval[0];
            //         std::cout<<std::endl<<"Level "<<level<<"All Cells "<<i;
            //     });

            //     projcells([&](auto, auto &interval, auto) {
            //         auto i = interval[0];
            //         std::cout<<std::endl<<"Level "<<level<<"Proj Cells "<<i;
            //     });
                
            //     cellsandghosts([&](auto, auto &interval, auto) {
            //         auto i = interval[0];
            //         std::cout<<std::endl<<"Level "<<level<<"Cells and ghosts "<<i;
            //     });
            // }


            // std::cout<<std::endl
            //          <<std::endl
            //          <<std::endl
            //          <<std::endl;

            // Initialization
            auto f  = init_f(mesh , 0.0);
            auto fR = init_f(meshR, 0.0);             

            double T = 0.4;
            double dx = 1.0 / (1 << max_level);
            double dt = dx;

            std::size_t N = static_cast<std::size_t>(T / dt);

            double t = 0.0;

            std::ofstream out_time_frames;
            std::ofstream out_error_exact_ref;
            std::ofstream out_diff_ref_adap;
            std::ofstream out_compression;

            out_time_frames.open     ("./d1q2/time_frame_s_"     +std::to_string(s)+"_eps_"+std::to_string(eps)+".dat");
            out_error_exact_ref.open ("./d1q2/error_exact_ref_"  +std::to_string(s)+"_eps_"+std::to_string(eps)+".dat");
            out_diff_ref_adap.open   ("./d1q2/diff_ref_adap_s_"  +std::to_string(s)+"_eps_"+std::to_string(eps)+".dat");
            out_compression.open     ("./d1q2/compression_s_"    +std::to_string(s)+"_eps_"+std::to_string(eps)+".dat");


            for (std::size_t nb_ite = 0; nb_ite < N; ++nb_ite)
            {
                tic();
                for (std::size_t i=0; i<max_level-min_level; ++i)
                {
                    //std::cout<<std::endl<<"Passe "<<i;
                    if (coarsening(f, eps, i))
                        break;
                }
                auto duration_coarsening = toc();

                // save_solution(f, eps, nb_ite, "coarsening");

                tic();
                for (std::size_t i=0; i<max_level-min_level; ++i)
                {
                    if (refinement(f, eps, 300.0, i))
                        break;
                }
                auto duration_refinement = toc();
                save_solution(f, eps, nb_ite, "refinement");


                //std::cout<<std::endl<<"Mesh before computing solution"<<std::endl<<mesh;


                // Create and initialize field containing the leaves
                tic();
                mure::Field<Config, int, 1> tag_leaf{"tag_leaf", mesh};
                tag_leaf.array().fill(0);
                mesh.for_each_cell([&](auto &cell) {
                    tag_leaf[cell] = static_cast<int>(1);
                });
                auto duration_leaf_checking = toc();

                mure::Field<Config, int, 1> tag_leafR{"tag_leafR", meshR};
                tag_leafR.array().fill(0);
                meshR.for_each_cell([&](auto &cell) {
                    tag_leafR[cell] = static_cast<int>(1);
                });

                auto error = compute_error(f, fR, t);

                out_time_frames    <<t       <<std::endl;
                out_error_exact_ref<<error[0]<<std::endl;
                out_diff_ref_adap  <<error[1]<<std::endl;
                out_compression    <<static_cast<double>(mesh.nb_cells(mure::MeshType::cells)) 
                                   / static_cast<double>(meshR.nb_cells(mure::MeshType::cells))<<std::endl;

                save_refined_solution(fR, min_level, max_level, eps, nb_ite);

                std::cout<<std::endl;

                

                
                tic();
                one_time_step(f, tag_leaf, s);
                auto duration_scheme = toc();

                tic();
                one_time_step(fR, tag_leafR, s);
                auto duration_schemeR = toc();

                t += dt;

                tic();
                save_solution(f, eps, nb_ite, "onetimestep");
                auto duration_save = toc();


                std::cout<<std::endl<<"\n=======Iteration "<<nb_ite<<"  time "<<t<<" summary========"
                                    <<"\nCoarsening: "<<duration_coarsening
                                    <<"\nRefinement: "<<duration_refinement
                                    <<"\nLeafChecking: "<<duration_leaf_checking
                                    <<"\nScheme: "<<duration_scheme
                                    <<"\nScheme reference: "<<duration_schemeR
                                    <<"\nSave: "<<duration_save
                                    <<"\nError exact - referece = "<< error[0]
                                    <<"\nError adaptive - referece = "<< error[1] << "\n";

                                    

            }
            
            out_time_frames.close();
            out_error_exact_ref.close();
            out_diff_ref_adap.close();
            out_compression.close();

        }

    }
    catch (const cxxopts::OptionException &e)
    {
        std::cout << options.help() << "\n";
    }



    return 0;
}
