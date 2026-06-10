#version 3.7;

#include "rad_def.inc"
#include "textures.inc"
#include "colors.inc"
#include "glass.inc"
#include "stones.inc"

// Right-handed coordinate system in which the z-axis points upwards
camera {
    // Aneurysm
    location <120,1200,120>
    look_at <45,0,45>
    // Wye
    location <100,240,225>
    look_at <25,20,45>
    sky z
    right -0.24*x*image_width/image_height
    up 0.24*z
}

global_settings
{
    assumed_gamma 1.0
    radiosity
    {
        Rad_Settings(Radiosity_Normal,off,off)
    }
    photons
    {
        count 20000
        media 200
        autostop 0
        jitter .4
    }
}

// White background
background{rgb 0.6}

// Two lights with slightly different colors
//light_source{<-80,-200,300> color rgb <0.77,0.75,0.75>}
//light_source{<350,-12,120> color rgb <0.38,0.40,0.40>}

light_source{<1000,1000,1000> color rgb <1.0,1.0,1.0>}
sphere
{
    <1000, 1000, 1000> 10
    texture
    {
        Glass
    }
    interior
    {
        ior 1.3
    }
    photons
    {
        target
        reflection on
        refraction on
    }
}

// Radius of the Voronoi cell network
#declare r=0.50;

// Radius of the particles
#declare s=1.00;

plane
{
    <0, 0, 1>, 0
    texture
    {
        pigment
        {
            Gray
        }
        finish
        {
            reflection 0.02
        }
    }
}

// Particles
union{
    #include "voronoi/generator_points_0000001.pov"
    texture
    {
        Ruby_Glass
    }
    interior
    {
        ior 1.3
    }
    photons
    {
        target
        reflection on
        refraction on
    }
}

// Voronoi cells
union{
    #include "voronoi/voronoi_cells_0000001.pov"
    texture
    {
        T_Stone13 
    }
    photons
    {
        target
        reflection on
    }
}
