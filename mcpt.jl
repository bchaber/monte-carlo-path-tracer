import LinearAlgebra: norm, cross, dot
Tuple3 = NTuple{3, Float64}
struct Ray
  origin :: Tuple3
  direction :: Tuple3
end
abstract type Material end
struct Diffusive  <: Material end
struct Specular   <: Material end
struct Refractive <: Material 
  εr :: Float64
end

struct Sphere
  radius :: Float64
  center :: Tuple3
  emission :: Tuple3
  color :: Tuple3
  material :: Material
end
corrected(x::Number) = floor(UInt8, (255.0clamp(x)^(1.0/2.2) + .5))
normalize(v::Tuple3) = v./sqrt(dot(v,v))
clamp(x::Number) = x < zero(x) ? zero(x) : x > one(x) ? one(x) : x
cross(u::Tuple3, v::Tuple3) = (u[2]*v[3]-u[3]*v[2], u[3]*v[1]-u[1]*v[3],u[1]*v[2]-u[2]*v[1])
dot(u::Tuple3,v::Tuple3) = sum(u .* v)

function intersect(ray::Ray, sph::Sphere)
  op = sph.center .- ray.origin
  b  = dot(op, ray.direction)
  det = b*b - dot(op, op) + sph.radius*sph.radius
  if det < 0. return 0.0 end
  det = sqrt(det)
  if (b - det) > 1e-4 return b - det end
  if (b + det) > 1e-4 return b + det end
  return 0.
end
const ex  = (1., 0., 0.)
const ey  = (0., 1., 0.)
const ez  = (0., 0., 1.)
const air = Refractive(1.0)
const spheres = (
  Sphere(1e5,  (+1e5+01,40.8,81.6),(0.,0.,0.), (.75, .25, .25), Diffusive()),
  Sphere(1e5,  (-1e5+99,40.8,81.6),(0.,0.,0.), (.25, .25, .75), Diffusive()),
  Sphere(1e5,  (50,40.8, 1e5),     (0.,0.,0.), (.75, .75, .75), Diffusive()),
  Sphere(1e5,  (50,40.8,-1e5+170), (0.,0.,0.), (.00, .00, .00), Diffusive()),
  Sphere(1e5,  (50, 1e5, 81.6),    (0.,0.,0.), (.75, .75, .75), Diffusive()),
  Sphere(1e5,  (50,-1e5+81.6,81.6),(0.,0.,0.), (.75, .75, .75), Diffusive()),
  Sphere(16.5, (27,16.5,47),       (0.,0.,0.), (.99, .99, .99), Specular()),
  Sphere(16.5, (73,16.5,78),       (0.,0.,0.), (.99, .99, .99), Refractive(1.5)),
  Sphere(600., (50,681.6-.27,81.6),(12.,12.,12.), (.00, .00, .00), Diffusive())
)

function intersect(ray::Ray)
  t = 1e20
  hit = nothing
  for sph=spheres
    d = intersect(ray, sph)
    if d > 0. && d < t
      t = d
      hit = sph
    end
  end
  return t, hit
end

function radiance(ray::Ray, depth::Integer)
  t, hit = intersect(ray)
  if hit == nothing
    return (0.,0.,0.)
  end
  x  = ray.origin .+ ray.direction .* t
  n  = normalize(x .- hit.center)
  nl = dot(n, ray.direction) < 0.0 ? n : -1.0.*n
  f  = hit.color
  p  = maximum(f)

  if depth > 4
    if rand() < p
      f = f./p
    end
    return hit.emission
  end
  return hit.emission .+ f.*radiance(hit.material, depth + 1, x, ray, n, nl)
end

function radiance(mat::Diffusive, depth, x, r, n, nl)
  r1 = rand()*2π
  r2 = rand()
  r2s= sqrt(r2)
  w = nl
  u = normalize(cross(abs(w[1]) > 0.1 ? ey : ex, w))
  v = cross(w, u)
  d = normalize(u.*cos(r1).*r2s .+ v.*sin(r1).*r2s .+ w.*sqrt(1.0-r2))

  return radiance(Ray(x,d), depth)
end

function radiance(mat::Specular, depth, x, r, n, nl)
  reflected = Ray(x, r.direction .- 2.0 .* dot(n, r.direction) .* n)
  return radiance(reflected, depth)
end

function radiance(mat::Refractive, depth, x, r, n, nl)
  reflected = Ray(x, r.direction .- 2.0 .* dot(n, r.direction) .* n) # Ideal dielectric REFRACTION
  into = dot(n, nl) > 0.0 # Ray from outside going in?
  nc = air.εr
  nt = mat.εr
  nnt= into ? nc/nt : nt/nc
  ddn= dot(r.direction, nl)
  cos2t = 1.0 - nnt*nnt*(1.0 - ddn*ddn)
  if cos2t < 0.0 # Total internal reflection
    return radiance(reflected, depth)
  end
  tdir = normalize(r.direction.*nnt .- (into ? n : -1.0.*n).*(ddn*nnt + sqrt(cos2t)))
  someray = Ray(x, tdir)
  a = nt - nc
  b = nt + nc
  R0= (a*a)/(b*b)
  c = 1.0 - (into ? -ddn : dot(tdir, n))
  Re= R0 + (1.0 - R0)*c*c*c*c*c
  Tr= 1.0 - Re
  P = .25 + .5Re
  RP= Re/P
  TP= Tr/(1.0-P)
  if depth > 2
    if rand() < P # Russian roulette
      return radiance(reflected, depth).*RP
    else
      return radiance(someray, depth).*TP
    end
  end
  return radiance(reflected, depth).*Re .+ radiance(someray, depth).*Tr
end

import Printf: @printf
function main(w=320, h=240, spp=1)
  cam = Ray((50.0, 52.0, 295.6), normalize((0.0, -0.042612, -1.0)))
  cx  = (.5135w/h, 0.0, 0.0)
  cy  = 0.5135.*normalize(cross(cx, cam.direction))
  c   = zeros(w, h, 3)
  d   = (0.0, 0.0, 0.0)
  r   = (0.0, 0.0, 0.0)

  for y=1:h # Loop over image rows
    @printf "\rRendering %d %5.2f%% " 4spp 100.0y/h
    for x=1:w # Loop cols
      for sy=1:2 # 2x2 subpixel rows
        r = (0.0, 0.0, 0.0)
        for sx=1:2 # 2x2 subpixel cols
          for s=1:spp
            r1 = 2.0rand()
            r2 = 2.0rand()
            dx = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1)
            dy = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2)
            d  = @. cx.*( ( (sx + .5 + dx)/2.0 + x)/w - .5) .+
                    cy.*( ( (sy + .5 + dy)/2.0 + y)/h - .5) .+
                    cam.direction;
            r = r .+ radiance(Ray(cam.origin .+ 140.0.*d, normalize(d)), 0)./spp
          end # Camera rays are pushed ^^^^^ forward to start in interior
          c[x,y,:].+= @. .25clamp(r)
        end
      end
    end
  end
  # Save image
  open("image.ppm", "w") do f
    write(f, "P6 $w $h 255\n")
    for y=reverse(1:h)
      for x=reverse(1:w)
        write(f, corrected.(c[x,y,:]))
      end
    end
  end
end
@time main(2,2,1)
import Profile: clear_malloc_data
clear_malloc_data()
#@time main(32,24,1)
@time main(320,240,125)
