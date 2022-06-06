#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl3.h>
//#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <string>
#include <iostream>
#include < cmath >
#define  SCREEN_WIDTH 1200
#define  SCREEN_HEIGHT 800
GLfloat V_user = 200;
GLfloat Feps_user = -100;
GLfloat D_user = 1;
GLfloat p_size = 10;
GLfloat a_user = 1;
GLfloat eps_user = p_size;
GLfloat rmax =5*p_size;
GLfloat delimetr =15;
GLfloat dt = 0;
GLfloat Turbulent_coeff = 0;

GLfloat linesize = 5;

GLfloat barrierx1 = 0;
GLfloat barriery1 = 0;
GLfloat barrierx2 = 0;
GLfloat barriery2 = 0;
GLfloat menysize = 250;

bool colormap_flag = true;


static void cursorPositionCallback(GLFWwindow* window, double xPos, double yPos);
static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);

bool stop_flag = false;
bool first_click = false;

class point {
public:
    GLfloat x=0;
    GLfloat y=0;
    GLfloat vx=0;
    GLfloat vy= 0;
    GLfloat m=0.1;
    bool was_intersec = false;
};

class barrier {
public:
    GLfloat x1 ;
    GLfloat y1 ;
    GLfloat x2 ;
    GLfloat y2 ;
};

std::vector<barrier> barriers;
std::vector<point> points;

GLvoid print_point(point point) {
    GLfloat pointVertex[] = { point.x,point.y };
    if (colormap_flag) {
        GLfloat current_v = sqrt(pow(point.vx, 2) + pow(point.vy, 2));
        GLfloat max_v = V_user * 2;
        GLfloat norm_v = current_v / max_v;
        glColor3f((norm_v), 0, (1 - norm_v));
    }
    else {
        glColor3f(1, 1, 1);
    }

    glEnable(GL_POINT_SMOOTH);
    glEnableClientState(GL_VERTEX_ARRAY);
    

    glPointSize(p_size);
    glVertexPointer(2, GL_FLOAT, 0, pointVertex);
    //glColorPointer(2, GL_FLOAT, 0, colour);
    glDrawArrays(GL_POINTS, 0, 1);
    glDisableClientState(GL_VERTEX_ARRAY);

    glDisable(GL_POINT_SMOOTH);
}

GLvoid print_line(barrier barrier) {
    GLfloat linetVertex[] = { barrier.x1,barrier.y1,barrier.x2,barrier.y2 };
    /*glEnable(GL_POINT_SMOOTH);*/
    glColor3f(255, 255, 255);
    glEnableClientState(GL_VERTEX_ARRAY);
    //glPointSize(p_size);
    glLineWidth(linesize);
    glVertexPointer(2, GL_FLOAT, 0, linetVertex);
    glDrawArrays(GL_LINES, 0, 2);
    glDisableClientState(GL_VERTEX_ARRAY);
    //glDisable(GL_POINT_SMOOTH);
}


GLvoid make_column(GLint num_points) {
    for (GLint i = GLint(1); i < num_points; i = i + GLint(1)) {
        points.push_back(point());
        points.back().x = SCREEN_WIDTH;
        points.back().y = SCREEN_HEIGHT / num_points * i + delimetr;
        points.back().vx = -V_user;
        points.back().vy = 2*(((double)rand() / (RAND_MAX)) - 0.5) * V_user * Turbulent_coeff;
    }
}


GLvoid Jones(int i, int j) {
   
    GLfloat rx = points[j].x - points[i].x;
    GLfloat ry = points[j].y - points[i].y;
    GLfloat r = sqrt(pow(rx,2.0)+ pow(ry, 2.0));
    //std::cout <<"rx ="<< rx << std::endl;
    //std::cout << "ry =" << ry << std::endl;
    //std::cout << "r =" << r << std::endl;
  
    GLfloat F = 0;
    if (r < eps_user) {
        F = Feps_user;
    }
    else {
        if (r > rmax) {
            F = 0;
        }
        else {
            F = (12 * D_user / a_user) * (pow(a_user / r, 7) - pow(a_user / r, 13));
        }
        
    }
    //std::cout << "F =" << F << std::endl;
    points[i].vx = points[i].vx + rx / r * F * dt / points[i].m;
    points[i].vy = points[i].vy + ry / r * F * dt / points[i].m;
    //points[i].vx = points[i].vx + rx / r * F * dt / points[i].m;
    //points[i].vy = points[i].vx + ry / r * F * dt / points[i].m;
    //std::cout << "vx =" << points[i].vx << std::endl;
    //std::cout << "vy =" << points[i].vy << std::endl;
}


bool ifintersec(point& point, barrier& barrier) {

    GLfloat localx1 = barrier.x1;
    GLfloat localx2 = barrier.x2;

    GLfloat localpointx = point.x ;
    GLfloat localpointy = point.y ;

    if (barrier.x2 < barrier.x1){
        localx1 = barrier.x2;
        localx2 = barrier.x1;


    }


    if (localx1 == localx2) {

        GLfloat localy1 = barrier.y1;
        GLfloat localy2 = barrier.y2;
        if (barrier.y2 < barrier.y1) {
            localy1 = barrier.y2;
            localy2 = barrier.y1;
        }


        if ((abs(localpointx + p_size / 2 - localx1) <= p_size/10)|| (abs(localpointx - p_size / 2 - localx1) <= p_size / 10)) {
            if ((localpointy - p_size / 2 < localy2) && (localpointy + p_size / 2 > localy1)) {
                return true;
            }
        }

        return false;
    }



    GLfloat cline = (barrier.y2- barrier.y1)/((barrier.x2 - barrier.x1));
    GLfloat climb = barrier.y1- cline* barrier.x1;



    GLfloat a = 1+pow(cline,2);
    GLfloat b = -2* point.x+2*cline*climb-2*cline* point.y;
    GLfloat c = pow(point.x, 2) + pow(climb, 2) + pow(point.y, 2) -2*climb* point.y - pow(p_size/2, 2);
    GLfloat Disc = pow(b, 2) - 4 * a * c;
    if (Disc < 0) {
        return false;
    }
    else {
        GLfloat ans1 = (-b + sqrt(Disc)) / (2 * a);
        GLfloat ans2 = (-b - sqrt(Disc)) / (2 * a);
        if ((localx1 <= ans1 && localx2 >= ans1)|| (localx1 <= ans2 && localx2 >= ans2)) {
            return true;
        }
            else {
                return false;
            }
        }


 }


GLvoid BarrierReflect(point& point, barrier& barrier) {
    GLfloat localx1 = barrier.x1;
    GLfloat localx2 = barrier.x2;
    GLfloat localy1 = barrier.y1;
    GLfloat localy2 = barrier.y2;
    if (barrier.x2 < barrier.x1) {
        localx1 = barrier.x2;
        localx2 = barrier.x1;
        localy1 = barrier.y2;
        localy2 = barrier.y1;
    }

    GLfloat cline = (barrier.y2 - barrier.y1) / ((barrier.x2 - barrier.x1));
    GLfloat climb = barrier.y1 - cline * barrier.x1;

    if (point.y == point.x*cline+climb) {
        point.vx = -point.vx;
        point.vy = -point.vy;
        return;
    }

    GLfloat len = sqrt(pow(localx2 - localx1, 2) + pow(localy2 - localy1, 2));
    GLfloat n_x = -(localy2 - localy1) / len;
    GLfloat n_y = (localx2 - localx1) / len;
    GLfloat projection = point.vx* n_x + point.vy * n_y;

    GLfloat reflected_vx = point.vx - 2 * n_x * projection;
    GLfloat reflected_vy = point.vy - 2 * n_y * projection;

    point.vx = reflected_vx;
    point.vy = reflected_vy;
}



GLvoid phys(GLint num_points) {
    
    if (points.back().x <= SCREEN_WIDTH - delimetr - p_size-30) {
        make_column(num_points);
    }



    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < points.size(); j++) {
            if (i != j) {
                Jones( i, j);
            }
        }
    }

    for (int i = 0; i < points.size(); i++) {
        //reflextion from bounds
        if (points[i].y + points[i].vy * dt-p_size<0 || points[i].y + points[i].vy * dt + p_size > SCREEN_HEIGHT) {
            points[i].vy = -points[i].vy;
        }
        //MOVE
        if ((points[i].x + points[i].vx * dt < 0 + p_size) || ((points[i].x + points[i].vx * dt > SCREEN_WIDTH + p_size))) {
            points.erase(points.begin() + i);
        }
        else {
            points[i].x = points[i].x + points[i].vx * dt;
            points[i].y = points[i].y + points[i].vy * dt;
        }

    }

    //reflextion from barriers
    
    for (int i = 0; i < points.size(); i++) {
        for (int j = 0; j < barriers.size(); j++) {
            //reflextion from bounds
            bool res_func_intersec = ifintersec(points[i], barriers[j]);
            if (res_func_intersec && points[i].was_intersec==false) {
                BarrierReflect(points[i], barriers[j]);
                points[i].was_intersec = true;
            }
            else if (res_func_intersec==false && points[i].was_intersec == true) {
                points[i].was_intersec = false;
                }
            
        }
    }

}




int main(void)
 
{
    /* Initialize the library */
    if (!glfwInit())
        return -1;

    std::cout << "Sucsessfully started!" << std::endl;
    std::cout << "Ver. 1.2" << std::endl;
    



    //num of points
    GLint num_points = int(std::floor(SCREEN_HEIGHT/(p_size+2*delimetr)));
    //create vec of points
    
    

    barriers.push_back(barrier());
    barriers[0].x1 = SCREEN_WIDTH*2/3;
    barriers[0].y1 = SCREEN_HEIGHT;
    barriers[0].x2 = SCREEN_WIDTH/2;
    barriers[0].y2 = SCREEN_HEIGHT-100;


    barriers.push_back(barrier());
    barriers[1].x1 = SCREEN_WIDTH * 2 / 3;
    barriers[1].y1 = 0;
    barriers[1].x2 = SCREEN_WIDTH / 2;
    barriers[1].y2 = 100;

    GLFWwindow* window;

    //create first column
    make_column(num_points);




    point mypoint;
    mypoint.x = 250;
    mypoint.y = 280;
    bool result = ifintersec( mypoint, barriers[0]);


    
    /* Create a windowed mode window and its OpenGL context */
    glfwWindowHint(GLFW_RESIZABLE, GLFW_FALSE);
    window = glfwCreateWindow(SCREEN_WIDTH+ menysize, SCREEN_HEIGHT, "Particles Flow", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }
    

    glfwSetCursorPosCallback(window, cursorPositionCallback );
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    glViewport(0.0f, 0.0f, SCREEN_WIDTH, SCREEN_HEIGHT); // specifies the part of the window to which OpenGL will draw (in pixels), convert from normalised to pixels
    glMatrixMode(GL_PROJECTION); // projection matrix defines the properties of the camera that views the objects in the world coordinate frame. Here you typically set the zoom factor, aspect ratio and the near and far clipping planes
    glLoadIdentity(); // replace the current matrix with the identity matrix and starts us a fresh because matrix transforms such as glOrpho and glRotate cumulate, basically puts us at (0, 0, 0)
    glOrtho(0, SCREEN_WIDTH, 0, SCREEN_HEIGHT, 0, 1); // essentially set coordinate system
    glMatrixMode(GL_MODELVIEW); // (default matrix mode) modelview matrix defines how your objects are transformed (meaning translation, rotation and scaling) in your world
    glLoadIdentity(); // same as above comment


        IMGUI_CHECKVERSION();

        ImGui::CreateContext();
     

        ImGuiIO& io = ImGui::GetIO(); (void)io;
        ImGui::StyleColorsDark();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330");


      

    /* Loop until the user closes the window */

    double last_time = 0;
    glfwSetTime(0);
    while (!glfwWindowShouldClose(window))
    {   

        num_points = int(std::floor(SCREEN_HEIGHT / (p_size + 2 * delimetr)));

      

        if (stop_flag) {
            double now = glfwGetTime();
            dt = now - last_time;
            last_time = now;
        }
        else {
            dt = 0;
        }
        glClearColor(0.07f, 0.13f, 0.17f, 1.0f);
        /*std::cout << "dt =" << dt << std::endl;*/
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);


        //GUI
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();


        phys(num_points);

        GLfloat pointVertex[] = { SCREEN_WIDTH/2,SCREEN_WIDTH/2 };

        //print_point(mypoint);
        /*for (std::vector<point>::iterator one_point = points.begin(); one_point != points.end(); one_point++) {*/
        for (int i = 0; i < points.size();i++) {
            print_point(points[i]);
        }
        for (int j = 0; j < barriers.size(); j++) {
            print_line(barriers[j]);
        }


        //GUI        
        ImGui::SetNextWindowPos(ImVec2(SCREEN_WIDTH,0));
        ImGui::SetNextWindowSize((ImVec2(menysize, SCREEN_HEIGHT)));

        ImVec2 scrollsize = ImVec2(10, 50);

        ImGui::Begin("Menu",0, ImGuiWindowFlags_NoMove);
        ImGui::Text("Mechanical:");

        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Size", &p_size,2.0f,30.0f);
        ImGui::SliderFloat("Delimetr", &delimetr, 5.0f, 30.0f);
        ImGui::SliderFloat("Velocity", &V_user, 1.0f, 500.0f);
        ImGui::SliderFloat("Turbulent", &Turbulent_coeff, 0.0f, 1.0f);
        ImGui::Text("Lennard-Jones parametres:");
        ImGui::SliderFloat("D", &D_user, 0.1f, 1.0f);
        ImGui::SliderFloat("a", &a_user, 0.1f, 50.0f);
        ImGui::SliderFloat("F min", &Feps_user, -100.0f, 0.0f);
        ImGui::SliderFloat("radius of interact", &rmax, 0.0f, 100.0f);
        
 

  
        ImGui::Text("Deleting:");
        if (ImGui::Button("Delete Particles")) {
            points = {};
           
        }

        if (ImGui::Button("Delete Barriers")) {
            barriers = {};
        }

        ImGui::Text("Control:");
       if (ImGui::Button("Stop")) {
            stop_flag=false;
        }


        if (ImGui::Button("Start")) {
            stop_flag = true;
            last_time = glfwGetTime();
        }
        ImGui::SetCursorPos(ImVec2(10, SCREEN_HEIGHT-30));
        ImGui::Checkbox("Colormap", &colormap_flag);


       

        //ImGui::SliderFloat("Size", &p_size, 1.0f, 30.0f);
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());


        /* Swap front and back buffers */
        glfwSwapBuffers(window);
        glClear(GL_COLOR_BUFFER_BIT);

        /* Poll for and process events */
        glfwPollEvents();
    }

    ImGui_ImplGlfw_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwTerminate();
    return 0;
}
static void cursorPositionCallback(GLFWwindow* window, double xPos, double yPos) {

};

static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT && action==GLFW_PRESS) {
       
        double xpos, ypos;
        //getting cursor position
        glfwGetCursorPos(window, &xpos, &ypos);
        if (xpos < SCREEN_WIDTH) {
            if (first_click == false) {
                barrierx1 = xpos;
                barriery1 = SCREEN_HEIGHT - ypos;
                first_click = true;
            }
            else {
                barrierx2 = xpos;
                barriery2 = SCREEN_HEIGHT - ypos;
                barriers.push_back(barrier());
                barriers.back().x1 = barrierx1;
                barriers.back().y1 = barriery1;

                barriers.back().x2 = barrierx2;
                barriers.back().y2 = barriery2;
                first_click = false;
            }
        }
    }
};