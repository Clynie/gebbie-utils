% Helper class for creating a movie
%
% Author: John Gebbie
% Creation Date: Oct 23, 2012

classdef MovieMaker < handle
    
    properties
        fps = 30
        res = [640 480]
        qual = 75
        fmt = 'Motion JPEG AVI'
        file = 'movie.avi'
        tmp_dir_suffix = '_frames'
        frame_num
    end
    
    methods
        function [o] = MovieMaker()
        end
        
        function [dir_name] = get_dir_name(o)
            dir_name = [o.file o.tmp_dir_suffix];
        end
        
        function [file_name] = get_file_name(o,frame_num)
            fn = sprintf('%09d.png',frame_num);
            file_name = fullfile(o.get_dir_name(), fn);
        end
        
        function [frame_nums] = get_existing_frame_nums(o)
            files = dir(fullfile(o.get_dir_name(),'*.png'));
            N = length(files);
            frame_nums = nan(N,1);
            for n = 1:N
                frame_nums(n) = str2double(files(n).name(1:end-4));
            end
            frame_nums = sort(frame_nums);
        end
        
        function [frame_num] = init_next_frame_num(o)
            % sets the next frame number
            frame_nums = o.get_existing_frame_nums();
            N = length(frame_nums);
            if N > 0
                o.frame_num = find(frame_nums ~= (1:N)',1,'first');
            else
                o.frame_num = 1;
            end
            frame_num = o.frame_num;
        end
        
        function [] = prep_figure(o)
            figure(1); close(1); figure(1);
            set(gcf,'Units','Pixels');
            tmp_pos = get(gcf,'Position');
            set(gcf,'Position',[tmp_pos(1:2) o.res]);
        end
        
        function [] = write_frame(o)
            if ~exist(o.get_dir_name(),'dir')
                mkdir(o.get_dir_name());
            end
            drawnow(); % flush drawing to screen
            im = frame2im(getframe(gcf));
            imwrite(im, fullfile(o.get_dir_name(), ...
                sprintf('%09d.png',o.frame_num)));
        end
        
        function [] = frames_to_movie(o)
            writer = VideoWriter(o.file,o.fmt);
            writer.FrameRate = o.fps;
            writer.Quality = o.qual;
            writer.open();
            try
                frame_nums = o.get_existing_frame_nums();
                for n = 1:length(frame_nums)
                    im = imread(o.get_file_name(frame_nums(n)));
                    writer.writeVideo(im);
                end
            catch ex
                writer.close();
                rethrow(ex);
            end
            writer.close();
        end
        
    end
end

