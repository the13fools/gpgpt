This is very similar to the 2d version, should combine the two projects probably but just forking for now.  


ffmpeg -r 5 -pattern_type glob -i 'cur_file_iter_*.png' -c:v libx264 high_res.mp4
