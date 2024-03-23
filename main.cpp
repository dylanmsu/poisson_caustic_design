/*#include "stdio.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "src/caustic_design.h"

#include <thread>

#define cimg_use_png
//#define cimg_use_jpeg

#if defined(_WIN32) || defined(_WIN64)
    #include "CImg.h"
#else
    #include<X11/Xlib.h>
    #include "/home/dylan/caustic_engineering/CImg.h"
#endif

void image_to_grid(const cimg_library::CImg<unsigned char>& image, std::vector<std::vector<double>>& image_grid) {
    for (int i = 0; i < image.height(); ++i) {
        std::vector<double> row;
        for (int j = 0; j < image.width(); ++j) {
            double r = image(j, i, 0) / 255.0; // Normalize R component
            double g = image(j, i, 1) / 255.0; // Normalize G component
            double b = image(j, i, 2) / 255.0; // Normalize B component
            double gray = (0.299 * r) + (0.587 * g) + (0.114 * b); // Calculate grayscale value using luminosity method
            row.push_back(gray);
        }
        image_grid.push_back(row);
    }
}

void save_grid_as_image(const std::vector<std::vector<double>>& img, int resolution_x, int resolution_y, const std::string& filename) {
    // Create an empty CImg object with the specified resolution
    cimg_library::CImg<unsigned char> image(resolution_x, resolution_y);

    // Copy the grid data to the image
    for (int i = 0; i < resolution_y; ++i) {
        for (int j = 0; j < resolution_x; ++j) {
            image(j, i) = static_cast<unsigned char>(img[i][j] * 255); // Scale to [0, 255] and cast to unsigned char
        }
    }

    // Save the image as a PNG
    image.save(filename.c_str());
}

std::unordered_map<std::string, std::string> parse_arguments(int argc, char const *argv[]) {
    // Define a map to store the parsed arguments
    std::unordered_map<std::string, std::string> args;

    // Iterate through command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string key, value;

        // Check if argument starts with '--'
        if (arg.substr(0, 2) == "--") {
            // Split argument by '=' to separate key and value
            size_t pos = arg.find('=');
            if (pos != std::string::npos) {
                key = arg.substr(2, pos - 2);
                value = arg.substr(pos + 1);
            }
            else {
                key = arg.substr(2);
                value = ""; // No value provided
            }
        }
        // Check if argument starts with '-'
        else if (arg[0] == '-') {
            // The next argument is the value
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                key = arg.substr(1);
                value = argv[++i];
            }
            else {
                key = arg.substr(1);
                value = ""; // No value provided
            }
        }
        // Invalid argument format
        else {
            //std::cerr << "Invalid argument format: " << arg << std::endl;
            //return 1;
        }

        // Store key-value pair in the map
        args[key] = value;
    }

    return args;
}*/

/*int main(int argc, char const *argv[])
{
    // Parse user arguments
    std::unordered_map<std::string, std::string> args = parse_arguments(argc, argv);

    // Load image
    cimg_library::CImg<unsigned char> image(args["intput_png"].c_str());
    double aspect_ratio = (double)image.width() / (double)image.height();

    // Print parsed arguments
    //for (const auto& pair : args) {
    //    //std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    //    printf("Key: %s, Value: %s\r\n", pair.first.c_str(), pair.second.c_str());
    //}

    Caustic_design caustic_design;

    int mesh_resolution_x = atoi(args["res_w"].c_str());
    double mesh_width = std::stod(args["width"]);

    caustic_design.set_mesh_resolution(mesh_resolution_x, mesh_resolution_x / aspect_ratio);
    caustic_design.set_domain_resolution(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio);
    caustic_design.set_mesh_size(mesh_width, mesh_width / aspect_ratio);

    caustic_design.set_lens_focal_length(std::stod(args["focal_l"]));
    caustic_design.set_lens_thickness(std::stod(args["thickness"]));
    caustic_design.set_solver_max_threads(atoi(args["max_threads"].c_str()));

    image = image.resize(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio, -100, -100, 3); // Resize using linear interpolation

    // Convert image to grid
    std::vector<std::vector<double>> pixels;
    
    image_to_grid(image, pixels);

    caustic_design.initialize_solvers(pixels);

    caustic_design.export_paramererization_to_svg("../parameterization_0.svg", 0.5f);
    
    for (int itr=0; itr<100; itr++) {
        printf("starting iteration %i\r\n", itr);
        
        double step_size = caustic_design.perform_transport_iteration();

        caustic_design.export_paramererization_to_svg("../parameterization_" + std::to_string(itr + 1) + ".svg", 0.5f);
        //caustic_design.export_inverted_transport_map("../inverted.svg", 0.5f);

        printf("step_size = %f\r\n", step_size);

        // convergence treshold for the parameterization
        if (step_size < 0.005) break;
    }

    printf("\033[0;32mTransport map solver done! Starting height solver.\033[0m\r\n");

    for (int itr=0; itr<3; itr++) {
        caustic_design.perform_height_map_iteration(itr);
    }

    printf("Height solver done! Exporting as solidified obj\r\n");

    //caustic_design.save_solid_obj_target("../output.obj");
    caustic_design.save_solid_obj_source("../output.obj");

    return 0;
}*/

#include <napi.h>

#include "src/caustic_design.h"

Caustic_design caustic_design;

int mesh_resolution_x = 100;
int mesh_resolution_y = 100;

double aspect_ratio = 1.0f;

double mesh_width = 1.0f;

int add(int x, int y){
 return (x+y);
}

Napi::Number loadImage(const Napi::CallbackInfo& info) 
{   
    Napi::Env env = info.Env();
    
    Napi::Array b = info[0].As<Napi::Array>();
    int width = info[1].As<Napi::Number>();
    int height = info[2].As<Napi::Number>();
    
    std::vector<double> pixels;
    for(int i=0; i<b.Length(); i++)
    {
        Napi::Value v = b[i];
        if (v.IsNumber())
        {
            double value = (double)v.As<Napi::Number>();
            //printf("value = %f\r\n", value);
            pixels.push_back(value);
        }
    }

    std::vector<std::vector<double>> pixels_2d;
    for (int y = 0; y < height; y++) {
        std::vector<double> row;
        for (int x = 0; x < width; x++) {
            row.push_back(pixels[y * width + x]);
        }
        pixels_2d.push_back(row);
    }

    caustic_design.load_image(pixels_2d);
    caustic_design.initialize_solvers();

    Napi::Number returnValue = Napi::Number::New(env, pixels.size());

    return returnValue;
}

Napi::Array flatten2DArray(const std::vector<std::vector<double>>& array2D, const Napi::Env& env) {
    Napi::Array result = Napi::Array::New(env);
    size_t index = 0;

    for (const auto& row : array2D) {
        for (const auto& element : row) {
            result[index++] = Napi::Number::New(env, element);
        }
    }

    return result;
}

Napi::Value getErrorGrid(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();

    // Assuming you have a 2D array 'pixels_2d' that you want to flatten
    std::vector<std::vector<double>> pixels_2d = caustic_design.get_error_grid();

    // Flatten the 2D array into a 1D array
    Napi::Array flattenedArray = flatten2DArray(pixels_2d, env);

    return flattenedArray;
}

Napi::Number runTransportIteration(const Napi::CallbackInfo& info) {
    Napi::Env env = info.Env();
    
    double step_size = caustic_design.perform_transport_iteration();

    Napi::Number returnValue = Napi::Number::New(env, step_size);

    return returnValue;
}

Napi::Number addWrapped(const Napi::CallbackInfo& info) {
  Napi::Env env = info.Env();

  //check if arguments are integer only.
  if(info.Length()<2 || !info[0].IsNumber() || !info[1].IsNumber()){
      Napi::TypeError::New(env, "arg1::Number, arg2::Number expected").ThrowAsJavaScriptException();
  }

  //convert javascripts datatype to c++
  Napi::Number first = info[0].As<Napi::Number>();
  Napi::Number second = info[1].As<Napi::Number>();

  //run c++ function return value and return it in javascript
  Napi::Number returnValue = Napi::Number::New(env, add(first.Int32Value(), second.Int32Value()));
  
  return returnValue;
}

Napi::Object Init(Napi::Env env, Napi::Object exports) 
{
    caustic_design.set_mesh_resolution(mesh_resolution_x, mesh_resolution_x / aspect_ratio);
    caustic_design.set_domain_resolution(4*mesh_resolution_x, 4*mesh_resolution_x / aspect_ratio);
    caustic_design.set_mesh_size(mesh_width, mesh_width / aspect_ratio);

    caustic_design.set_lens_focal_length(1.0f);
    caustic_design.set_lens_thickness(0.2f);
    caustic_design.set_solver_max_threads(8);

    //export Napi function
    exports.Set("add", Napi::Function::New(env, addWrapped));
    exports.Set("loadImage", Napi::Function::New(env, loadImage));
    exports.Set("getErrorGrid", Napi::Function::New(env, getErrorGrid));
    exports.Set("runTransportIteration", Napi::Function::New(env, runTransportIteration));
    return exports;
}

NODE_API_MODULE(addon, Init)
