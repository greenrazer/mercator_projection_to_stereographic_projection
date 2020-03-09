import numpy as np
from PIL import Image
from tqdm import tqdm
import argparse

DEGREES_PER_RADIAN = 180/np.pi
RADIANS_PER_DEGREE = np.pi/180

def deg_to_rad(deg):
  return deg*RADIANS_PER_DEGREE

def rad_to_deg(rad):
  return rad*DEGREES_PER_RADIAN

def rad_to_x_y(rx,ry, width, height):

  dx = rad_to_deg(rx)
  dy = rad_to_deg(ry)

  # Make sure to bound dy to avoid drawing over Antarctica
  dy = np.clip(dy,-180,179)
  
  # pixels per radian for long and lat
  inc = width/360

  # To pixel value
  px = dx*inc + width/2
  py = dy*inc + height/2

  # Bound x and y if needed
  x = np.round(px) % (width-1)
  y = np.round(py) % (height-1)

  return int(x), int(y) 

# https://mathworld.wolfram.com/MercatorProjection.html
def merc_long_lat_to_x_y(long_lambda, lat_phi):

  # convert from degrees to radians
  long_rad = deg_to_rad(long_lambda) # -pi to pi
  lat_rad = deg_to_rad(lat_phi) # -pi/2 to pi/2

  # Mercator projection
  merc_x = long_rad
  sin_phi = np.sin(lat_rad)
  merc_y = np.log((1+sin_phi)/(1-sin_phi))/2

  return merc_x, merc_y

# https://mathworld.wolfram.com/StereographicProjection.html
def ster_x_y_to_long_lat(x, y, off_long, off_lat):

  off_long_rad = deg_to_rad(off_long)
  off_lat_rad = deg_to_rad(off_lat)

  R = 1
  rho = np.sqrt(np.square(x)+np.square(y))
  c = 2*np.arctan2(rho,2*R)

  sin_c = np.sin(c)
  cos_c = np.cos(c)
  sin_off_lat = np.sin(off_lat_rad)
  cos_off_lat = np.cos(off_lat_rad)

  long_num = x*sin_c
  long_denom = rho*cos_off_lat*cos_c - y*sin_off_lat*sin_c
  long_rad = off_long_rad + np.arctan2(long_num, long_denom) # -pi to pi

  lat_in_1 = cos_c*sin_off_lat
  lat_in_2_num = y*sin_c*cos_off_lat
  lat_in_2 = 0 if lat_in_2_num == 0 else lat_in_2_num/rho
  lat_rad = np.arcsin(lat_in_1 + lat_in_2)

  # convert from radians to degrees
  long_lambda = rad_to_deg(long_rad)
  lat_phi = rad_to_deg(lat_rad)

  return long_lambda, lat_phi

def main(input_file, out_width, out_height, target_lat, target_long, scale_factor, out_file):

  out = Image.new('RGB', (out_width, out_height))
  out_pixels = out.load()

  scale = np.sqrt(np.square(out_width)+np.square(out_height))*0.05*scale_factor


  with Image.open(input_file) as im:
    in_width, in_height = im.size
    merc_pixels = im.load()

    with tqdm(total=out_width*out_height) as pbar:
      for x in range(out_width):
          for y in range(out_height):

              #Scale and translate point
              tx = (x - out_width/2)/scale
              ty = (out_height/2 - y)/scale

              #Inverse Sterographic Projction
              lon, lat = ster_x_y_to_long_lat(tx, ty, target_long, target_lat)

              #Mercator projection
              merc_rx, merc_ry = merc_long_lat_to_x_y(lon,lat)

              #Radians to pixel values
              mx, my = rad_to_x_y(merc_rx, -merc_ry, in_width, in_height)

              out_pixels[x,y] = merc_pixels[mx, my]
              pbar.update(1)
  if out_file:
    out.save(out_file)
  else:
    out.show()

if __name__ == '__main__':

  lat_bounds = (-90,90)
  long_bounds = (-180,180)
  
  parser = argparse.ArgumentParser(description='Convert Mercador projection map into sterographic projection map.')
  parser.add_argument('input_filename', type=str, help='Mercator projection image filename.')
  parser.add_argument('-wi', '--width', type=int, default=100, help='Width of the output map.')
  parser.add_argument('-he', '--height', type=int, default=100, help='Height of the output map.')
  parser.add_argument('-lat', '--latitude', type=float, default=0, help=f'Center latitude. Between {lat_bounds[0]} and {lat_bounds[1]}(Suffix "N" is positive, suffix "S" is negitive).')
  parser.add_argument('-long', '--longitude', type=float, default=0, help=f'Center longitude. Between {long_bounds[0]} and {long_bounds[1]}(Suffix "E" is positive, suffix "W" is negitive).')
  parser.add_argument('-sf', '--scale_factor', type=float, default=1, help=f'Scale factor. Zoom of projection.')
  parser.add_argument('-f', '--out_file', type=str, default=None, help='Directory to put output images in.')
  args = parser.parse_args()

  if args.width < 0:
    raise argparse.ArgumentTypeError(f"Width must be positive.")

  if args.height < 0:
    raise argparse.ArgumentTypeError(f"Height must be positive.")

  if lat_bounds[0] > args.latitude or args.latitude > lat_bounds[1]:
    raise argparse.ArgumentTypeError(f"Latitude must be between {lat_bounds[0]} and {lat_bounds[1]}.")

  if long_bounds[0] > args.longitude or args.longitude > long_bounds[1]:
    raise argparse.ArgumentTypeError(f"Longitude must be between {long_bounds[0]} and {long_bounds[1]}.")

  if args.scale_factor < 0:
    raise argparse.ArgumentTypeError(f"Scale factor must be positive.")

  main(args.input_filename, args.width, args.height, args.latitude, args.longitude, args.scale_factor, args.out_file)