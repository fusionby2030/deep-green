#include "raylib.h"
#include <math.h>
#include <vector> 

// Volvo Tundra Concept
Color volvo_butterscotch = {217, 154, 77, 255}; 
Color volvo_black = {40, 42, 40, 255};
Color volvo_ivory = {248, 254, 235, 255}; 
Color volvo_green = {79, 179, 82, 255}; 
Color volvo_green_muted = {79, 179, 82, 75};

const std::vector<float> temperature_outside = {-6.199999999999989, -6.399999999999977, -6.699999999999989, -6.699999999999989, -6.699999999999989, -6.800000000000011, -6.800000000000011, -7.100000000000023, -7.0, -7.600000000000023, -7.699999999999989, -8.199999999999989, -8.699999999999989, -8.600000000000023, -8.600000000000023, -8.800000000000011, -9.0, -8.899999999999977, -8.800000000000011, -8.600000000000023, -8.899999999999977, -9.0, -9.100000000000023, -8.699999999999989, -9.100000000000023, -9.300000000000011, -9.600000000000023, -9.5, -9.100000000000023, -9.199999999999989, -9.199999999999989, -9.199999999999989, -9.300000000000011, -9.600000000000023, -9.5, -9.600000000000023, -9.800000000000011, -9.800000000000011, -9.800000000000011, -9.600000000000023, -9.800000000000011, -10.0, -9.699999999999989, -10.100000000000023, -10.0, -10.199999999999989, -9.899999999999977, -9.899999999999977, -10.0, -9.800000000000011, -9.699999999999989, -9.699999999999989, -9.800000000000011, -9.800000000000011, -9.5, -8.600000000000023, -7.600000000000023, -7.399999999999977, -6.800000000000011, -6.300000000000011, -5.899999999999977, -5.699999999999989, -5.199999999999989, -4.600000000000023, -4.0, -3.8000000000000114, -3.3999999999999773, -3.0, -2.8000000000000114, -2.6000000000000227, -2.3000000000000114, -2.0, -1.8000000000000114, -1.5, -1.3999999999999773, -1.1999999999999886, -1.1000000000000227, -1.0, -1.0, -0.8999999999999773, -0.8000000000000114, -0.6999999999999886, -0.5, -0.39999999999997726, -0.5, -0.10000000000002274, -0.10000000000002274, 0.39999999999997726, 0.30000000000001137, -0.10000000000002274, 0.19999999999998863, 0.19999999999998863, -0.10000000000002274, 0.10000000000002274, -0.10000000000002274, -0.19999999999998863, -0.10000000000002274, -0.19999999999998863, -0.5, -0.6999999999999886, -1.3000000000000114, -1.8999999999999773, -2.1000000000000227, -2.3999999999999773, -2.8000000000000114, -2.8999999999999773, -3.3000000000000114, -3.8000000000000114, -3.8999999999999773, -4.100000000000023, -4.100000000000023, -4.100000000000023, -4.199999999999989, -4.100000000000023, -4.399999999999977, -4.800000000000011, -4.699999999999989, -4.699999999999989, -5.100000000000023, -5.399999999999977, -5.100000000000023, -5.199999999999989, -5.199999999999989, -5.300000000000011, -5.199999999999989, -5.199999999999989, -5.5, -5.600000000000023, -5.5, -5.399999999999977, -5.699999999999989, -5.899999999999977, -5.899999999999977, -6.100000000000023, -5.600000000000023, -6.0, -6.600000000000023, -6.199999999999989, -6.199999999999989, -6.300000000000011, -6.600000000000023, -6.5, -6.300000000000011, -6.600000000000023};


// Create example vector of data 
/* Maris ideas
- Check boxes for which data to use
    - Sunlight
    - Temperature outside
    - Temperature inside
- Day vs Week vs Month vs Year check box

*/ 

void draw_time_series_background_old(Rectangle rec) {
    DrawRectangleLinesEx(rec, 2.0, volvo_green_muted);

    // rec.x is the starting point of the x-axis
    // rec.y is the starting point of the y-axis
    float real_start_time = 0.0; 
    float real_end_time = 60.0*60.0*24.0; // 24 hours in seconds

    // Draw Vertical lines
    int num_ver_lines = 12;
    for (int i = 0; i < num_ver_lines; i++) {
        DrawLine(rec.x + i * rec.width / num_ver_lines, rec.y, rec.x + i * rec.width / num_ver_lines, rec.y + rec.height, volvo_green_muted);
        // Label the vertical lines by hour (intervals of 2 hours)
        DrawText(TextFormat("%d", i*2), rec.x + i * rec.width / num_ver_lines - 5, rec.y + rec.height + 10, 30, volvo_ivory); 
    }
    DrawText(TextFormat("%d", 24), rec.x + 24 * rec.width / 24 - 10, rec.y + rec.height + 10, 30, volvo_ivory); 

    // Plot temperature outside from the vector 
    float llim_y = -10.0; 
    float ulim_y = 10.0; 
    float delta_t = 600.0; // 10 minutes in seconds, temperatur epoints are sampled with this rate 

    for (size_t i = 1; i < temperature_outside.size(); i++) {
        float x = rec.x + i * rec.width / temperature_outside.size(); 
        float y = rec.y + rec.height - (temperature_outside.at(i) - llim_y) / (ulim_y - llim_y) * rec.height; 
        DrawCircle(x, y, 2, volvo_ivory);
    }
    // Draw Horizontal lines
    // Should represent evenly spaced temperature values from llim_y to ulim_y 
    int num_hor_lines = 10;
    for (int i = 0; i < num_hor_lines; i++) {
        DrawLine(rec.x, rec.y + i * rec.height / num_hor_lines, rec.x + rec.width, rec.y + i * rec.height / num_hor_lines, volvo_green_muted);
        // Label the horizontal lines by temperature
        int temp_value = (int)ulim_y - i * (int)(ulim_y - llim_y) / num_hor_lines;
        DrawText(TextFormat("%3d", temp_value), rec.x - 60, rec.y + i * rec.height / num_hor_lines - 10, 30, volvo_butterscotch);
    }
    DrawText(TextFormat("%3d", (int)llim_y), rec.x - 60, rec.y + rec.height - 10, 30, volvo_butterscotch);
}


void draw_time_series_background(Rectangle rec, int time_frame){
    float llim_y = -10.0; 
    float ulim_y = 10.0; 

    // Draw Horizontal lines
    // Should represent evenly spaced temperature values from llim_y to ulim_y 
    int num_hor_lines = 10;
    for (int i = 0; i < num_hor_lines; i++) {
        DrawLine(rec.x, rec.y + i * rec.height / num_hor_lines, rec.x + rec.width, rec.y + i * rec.height / num_hor_lines, volvo_green_muted);
        // Label the horizontal lines by temperature
        int temp_value = (int)ulim_y - i * (int)(ulim_y - llim_y) / num_hor_lines;
        DrawText(TextFormat("%3d", temp_value), rec.x - 60, rec.y + i * rec.height / num_hor_lines - 10, 30, volvo_butterscotch);
    }
    DrawText(TextFormat("%3d", (int)llim_y), rec.x - 60, rec.y + rec.height - 10, 30, volvo_butterscotch);

    // Draw Vertical lines
    // Text below changes based on time_frame 
    // 0 -> 24 hours, 1 -> 1 week, 2 -> 1 month, 3 -> 1 year
    // 0 -> with 12 lines gives 2 hours intervals
    // 1 -> with 14 lines gives 0.5 day intervals
    // 2 -> with 15 lines gives ~2 days intervals
    // 3 -> with 12 lines gives ~30 days intervals
    int num_ver_lines;
    float scaling_factor;
    int max_label;
    switch (time_frame) {
        case 0: 
            num_ver_lines = 12;
            scaling_factor = 2; 
            max_label = 24;
            break;
        case 1:
            num_ver_lines = 7;
            scaling_factor = 1; 
            max_label = 7;
            break;
        case 2:
            num_ver_lines = 15;
            scaling_factor = 2;
            max_label = 30;
            break;
        case 3:
            num_ver_lines = 12;
            scaling_factor = 1;
            max_label = 12;
            break;
        default:
            num_ver_lines = 12;
            scaling_factor = 2;
            max_label = 24;
            break;
    }
    for (int i = 0; i < num_ver_lines; i++) {
        DrawLine(rec.x + i * rec.width / num_ver_lines, rec.y, rec.x + i * rec.width / num_ver_lines, rec.y + rec.height, volvo_green_muted);
        // Label vertical lines based on time_frame
        DrawText(TextFormat("%d", (int)(i * scaling_factor)), rec.x + i * rec.width / num_ver_lines - 5, rec.y + rec.height + 10, 30, volvo_ivory);
    }
    DrawText(TextFormat("%d", max_label), rec.x + rec.width - 10, rec.y + rec.height + 10, 30, volvo_ivory);


    DrawRectangleLinesEx(rec, 2.0, volvo_green_muted);
}

void draw_measured_data(Rectangle rec, std::vector<float> &meas, int time_frame) {
    // Draw measured data based on the selected time frame
    // time_frame: 0 -> 24 hours, 1 -> 1 week, 2 -> 1 month, 3 -> 1 year
    float llim_y = -10.0; 
    float ulim_y = 10.0; 
    float delta_t = 60.0; // 1 minute in seconds, measured data is once per minute 

    size_t num_points_per_day = 60 * 24; // 60 minutes per hour * 24 hours
    size_t num_points;

    switch (time_frame) {
        case 0: // 24 hours
            num_points = num_points_per_day;
            break;
        case 1: // 1 week
            num_points = num_points_per_day * 7;
            break;
        case 2: // 1 month (assuming 30 days)
            num_points = num_points_per_day * 30;
            break;
        case 3: // 1 year (assuming 365 days)
            num_points = num_points_per_day * 365;
            break;
        default:
            num_points = num_points_per_day;
            break;
    }

    size_t start_index = meas.size() > num_points ? meas.size() - num_points : 0;

    for (size_t i = start_index; i < meas.size(); i++) {
        float x = rec.x + (i - start_index) * rec.width / num_points; 
        float y = rec.y + rec.height - (meas.at(i) - llim_y) / (ulim_y - llim_y) * rec.height; 
        DrawCircle(x, y, 2, volvo_butterscotch);
    }
}


void load_measured_data(std::vector<float> &meas){
    // create fake data
    // sampled every second, so sine wave with period of 24 hours 
    for (size_t i = 0; i < meas.size(); i++) {
        meas.at(i) = 5.0 * sin(2.0 * 3.1415 * i / (60.0*24.0));
    }


}

int main(void)  {
    int screenWidth = 1600;
    int screenHeight = 800;
    char text[] = "Tervetuloa Luovaklubbi!";
    InitWindow(screenWidth, screenHeight, "GHOUSE!");
    Font font = LoadFontEx("/usr/share/fonts/NationalPark-VariableVF.ttf", 96, 0, 0);

    GenTextureMipmaps(&font.texture);
    float fontSize = (float)font.baseSize - 2.0f;
    Vector2 fontPosition = { 20.0f, screenHeight/2.0f - 80.0f };
    Vector2 textSize = { 0.0f, 0.0f };
    
    textSize = MeasureTextEx(font, text, fontSize, 0);
    SetTextureFilter(font.texture, TEXTURE_FILTER_BILINEAR);
    int currentFontFilter = 0;      // TEXTURE_FILTER_POINT

    SetTargetFPS(60);               // Set our game to run at 60 frames-per-second
    // Rectangle textrec = {fontPosition.x, fontPosition.y, textSize.x + 0.5f, textSize.y+ 0.5f};

    Rectangle main_timeseries_rec = {300, 100,screenWidth / 2, screenHeight - 200}; 

    // create empty vector for measured data, 60*24*365 elements 
    std::vector<float> measured_data(60*24*365, 0.0);
    load_measured_data(measured_data);

    int time_frame = 1; // 0 -> 24 hours, 1 -> 1 week, 2 -> 1 month, 3 -> 1 year
    bool redraw = true;
    const char* time_frame_str = "Week";
    
    while (!WindowShouldClose()) {
        
        BeginDrawing();
            
            // DrawTextEx(font, "DASHBOARD TEXT", {fontPosition.x, fontPosition.y + 100}, fontSize, 1.0, volvo_green_muted);
            // draw_time_series_background(main_timeseries_rec, time_frame);
            // DrawTextEx(font, "the FUTURE", {main_timeseries_rec.x, main_timeseries_rec.y - fontSize}, fontSize, 1.0, volvo_ivory);
            // draw_measured_data(main_timeseries_rec, measured_data, time_frame);

            if (IsKeyPressed(KEY_ONE)) {
                time_frame = 0;
                redraw = true;
                time_frame_str = "Day"; 
            }
            if (IsKeyPressed(KEY_TWO)) {
                time_frame = 1;
                redraw = true;
                time_frame_str = "Week"; 
            }
            if (IsKeyPressed(KEY_THREE)) {
                time_frame = 2;
                redraw = true;
                time_frame_str = "Month"; 
            }
            if (IsKeyPressed(KEY_FOUR)) {
                time_frame = 3;
                redraw = true;                
                time_frame_str = "Year"; 
            }
            if (redraw) {
                ClearBackground(volvo_black);
                draw_time_series_background(main_timeseries_rec, time_frame);
                draw_measured_data(main_timeseries_rec, measured_data, time_frame);
                redraw = false;
                DrawTextEx(font, time_frame_str, {main_timeseries_rec.x, main_timeseries_rec.y - fontSize}, fontSize, 1.0, volvo_ivory);
            }

        EndDrawing();
    }
    CloseWindow(); 

    return 0;
}
