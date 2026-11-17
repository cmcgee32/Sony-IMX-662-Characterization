function imx662_caption(fig, ~, ~)
%IMX662_CAPTION  Apply the academic figure style (no explanatory text box).
%   Kept with the old signature so existing calls still work; the second
%   and third arguments (old "How to read this" text) are ignored.
    figure(fig);
    imx662_style(fig);
end
