#include "Triangle.h"
#include <GL/glew.h>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>

// A function clamping the input values to the lower and higher bounds
#define CLAMP(in, low, high) ((in) < (low) ? (low) : ((in) > (high) ? (high) : in))

using std::sqrt;
using std::pow;
using std::min;
using std::max;
using std::log2;

Triangle::Triangle()
{
	v[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	v[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	c[0] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[1] = glm::vec3(0.0f, 0.0f, 0.0f);
	c[2] = glm::vec3(0.0f, 0.0f, 0.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);
}

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);

	t[0] = glm::vec2(0.0f, 0.0f);
	t[1] = glm::vec2(0.0f, 0.0f);
	t[2] = glm::vec2(0.0f, 0.0f);

};

Triangle::Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2)
{
	v[0] = v0;
	v[1] = v1;
	v[2] = v2;

	t[0] = t0;
	t[1] = t1;
	t[2] = t2;

	c[0] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[1] = glm::vec3(1.0f, 1.0f, 1.0f);
	c[2] = glm::vec3(1.0f, 1.0f, 1.0f);

};

void Triangle::setColor(int vert, glm::vec3& color){
	c[vert].x = color.x;
	c[vert].y = color.y;
	c[vert].z = color.z;
}

glm::vec3 Triangle::getVert(int vert){
	return v[vert];
}

glm::vec3 getBarycentricConstants(glm::vec3& pos, glm::vec4 (&v)[3]){
	float x = pos.x;
	float y = pos.y;
	float x_a = v[0].x; float y_a = v[0].y;
	float x_b = v[1].x; float y_b = v[1].y;
	float x_c = v[2].x; float y_c = v[2].y;

	float alpha = ((-(x - x_b) * (y_c - y_b)) + ((y - y_b) * (x_c - x_b))) / ((-(x_a - x_b) * (y_c - y_b)) + ((y_a - y_b) * (x_c - x_b)));
	float beta = ((-(x - x_c) * (y_a - y_c)) + ((y - y_c) * (x_a - x_c))) / ((-(x_b - x_c) * (y_a - y_c)) + ((y_b - y_c) * (x_a - x_c)));
	float gamma = 1 - alpha - beta;
	glm::vec3 res(alpha, beta, gamma);
	return res;
}

// Rendering the triangle using OpenGL
void Triangle::RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool isTextured)
{

	glMatrixMode(GL_MODELVIEW);
	glLoadMatrixf(glm::value_ptr(modelViewMatrix));

	glMatrixMode(GL_PROJECTION);
	glLoadMatrixf(glm::value_ptr(projectionMatrix));
	
	// For textured object
	if (isTextured)
	{
		glEnable(GL_TEXTURE_2D);

		// Avoid modulating the texture by vertex color
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

		glBegin(GL_TRIANGLES);

			glTexCoord2f(t[0].x, t[0].y);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glTexCoord2f(t[1].x, t[1].y);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glTexCoord2f(t[2].x, t[2].y);
			glVertex3f(v[2].x, v[2].y, v[2].z);

		glEnd();

		glDisable(GL_TEXTURE_2D);


	}
	// For object with only vertex color
	else
	{
		glBegin(GL_TRIANGLES);

			glColor3f(c[0].x, c[0].y, c[0].z);
			glVertex3f(v[0].x, v[0].y, v[0].z);

			glColor3f(c[1].x, c[1].y, c[1].z);
			glVertex3f(v[1].x, v[1].y, v[1].z);

			glColor3f(c[2].x, c[2].y, c[2].z);
			glVertex3f(v[2].x, v[2].y, v[2].z);
		
		glEnd();
	}

}

float wrap(float x, float mx){
	if(x < 0){
		x += mx;
	}
	if(x > mx){
		x -= mx;
	}
	return x;
}

void getTexCoords(glm::vec3 &curr_px, glm::vec4 (&ptsHomo)[3], glm::vec3 &z_inv, glm::vec2 (&Q_sca)[3], float &res_x, float &res_y){
	glm::vec3 bary = getBarycentricConstants(curr_px, ptsHomo);
	float z_inv_interp = bary.x * z_inv[0] + bary.y * z_inv[1] + bary.z * z_inv[2];
	glm::vec2 Q_sca_interp = bary.x * Q_sca[0] + bary.y * Q_sca[1] + bary.z * Q_sca[2];
	Q_sca_interp /= z_inv_interp;
	float tex_x = Q_sca_interp.x;
	float tex_y = Q_sca_interp.y;
	// wraparound boundary condition
	tex_x = wrap(tex_x, 1);
	tex_y = wrap(tex_y, 1);
	res_x = tex_x;
	res_y = tex_y;
}

// Render the triangle on CPU
void Triangle::RenderCPU(
	glm::mat4 &modelViewMatrix, 
	glm::mat4 &projectionMatrix, 
	float height, // window height
	float width,  // window width
	float *color,
	std::vector<float*> textures, // texture images
	std::vector<size_t> tex_height, // heights of texture images
	std::vector<size_t> tex_width,  // widths of texture images
	int tex_mode, // 0: NN, 1: bilinear interp, 2: mipmapping
	float *depth
){
		auto at = [&](float *ptr, size_t r, size_t c, int col = -1) {
			int n = height;
			int m = width;
			if(col == -1){
				return *(ptr + (r * m) + c);
			} else {
				return *(ptr + (r * (m * 3) + (c * 3) + col));
			}
		};
		auto at_tex = [&](float *ptr, size_t r, size_t c, int lvl, int col = -1){
			int n = tex_height[lvl];
			int m = tex_width[lvl];
			if(col == -1){
				return *(ptr + (r * m) + c);
			} else {
				return *(ptr + (r * (m * 3) + (c * 3) + col));
			}
		};
		// wraps around edges of texture image
		auto at_tex_wrap = [&](float *ptr, size_t r, size_t c, int lvl, int col = -1){
			int n = tex_height[0];
			int m = tex_width[0];
			// r = WRAP()
			if(r >= n){
				r = 0;
			}
			if(c >= m){
				c = 0;
			}
			if(r < 0){
				r = n - 1;
			}
			if(c < 0){
				c = m - 1;
			}
			if(col == -1){
				return *(ptr + (r * m) + c);
			} else {
				return *(ptr + (r * (m * 3) + (c * 3) + col));
			}
		};
		auto set = [&](float *ptr, int r, int c, int col = -1, float put = -1){
			int n = height;
			int m = width;
			if(col == -1){
				*(ptr + (r * m) + c) = put;
			} else {
				*(ptr + (r * (m * 3) + (c * 3) + col)) = put;
			}
		};

		auto bilinear_interp = [&](int D, size_t img_x_close, size_t img_y_close, float tex_x, float tex_y){
			// get neighboring texels
			glm::vec3 all[4] = {
				glm::vec3(
					at_tex_wrap(textures[D], img_y_close, img_x_close, D, 0), 
					at_tex_wrap(textures[D], img_y_close, img_x_close, D, 1), 
					at_tex_wrap(textures[D], img_y_close, img_x_close, D, 2)
				),
				glm::vec3(
					at_tex_wrap(textures[D], img_y_close, img_x_close - 1, D, 0),
					at_tex_wrap(textures[D], img_y_close, img_x_close - 1, D, 1),
					at_tex_wrap(textures[D], img_y_close, img_x_close - 1, D, 2)
				),
				glm::vec3(
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close, D, 0),
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close, D, 1),
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close, D, 2)
				),
				glm::vec3(
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close - 1, D, 0),
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close - 1, D, 1),
					at_tex_wrap(textures[D], img_y_close - 1, img_x_close - 1, D, 2)
				)
			};
			float img_x, img_y;
			img_x = tex_x * tex_width[0];
			img_y = tex_y * tex_height[0];
			// compute interp weights and interp
			float w_1 = (img_x - img_x_close + 1) / 2.0f;
			float w_2 = 1.0f - w_1;
			glm::vec3 y1 = all[0] * w_1 + all[1] * w_2; 
			glm::vec3 y2 = all[2] * w_1 + all[3] * w_2;
			w_1 = (img_y - img_y_close + 1) / 2.0f;
			w_2 = 1.0f - w_1;
			glm::vec3 res = y1 * w_1 + y2 * w_2;
			return res;
		};

		glm::vec4 ptsHomo[3] = {
			glm::vec4(v[0].x, v[0].y, v[0].z, 1), 
			glm::vec4(v[1].x, v[1].y, v[1].z, 1), 
			glm::vec4(v[2].x, v[2].y, v[2].z, 1)
		};
		// apply modelview + projection transform to each vertex
		for(int i = 0; i < 3; i++){
			ptsHomo[i] = modelViewMatrix * ptsHomo[i];
			ptsHomo[i] = projectionMatrix * ptsHomo[i];
			ptsHomo[i].x /= ptsHomo[i].w;
			ptsHomo[i].y /= ptsHomo[i].w;
			// ptsHomo[i].z /= ptsHomo[i].w;
			ptsHomo[i].w /= ptsHomo[i].w;
		}
		// z_inv and Q_inv used to perspective correct tex coord interp
		glm::vec3 z_inv = {
			1.0f / ptsHomo[0].z,
			1.0f / ptsHomo[1].z,
			1.0f / ptsHomo[2].z
		};
		
		glm::vec2 Q_sca[3] = {
			glm::vec2(t[0].x / ptsHomo[0].z, t[0].y / ptsHomo[0].z),
			glm::vec2(t[1].x / ptsHomo[1].z, t[1].y / ptsHomo[1].z),
			glm::vec2(t[2].x / ptsHomo[2].z, t[2].y / ptsHomo[2].z)
		};
		// bounding box
		float minX, minY;
		float maxX, maxY;
		minX = minY = 1e9;
		maxX = maxY = -1e9;
		for(int i = 0; i < 3; i++){
			// std::cout << ptsHomo[i].x << " " << ptsHomo[i].y << " " << ptsHomo[i].z << " " << ptsHomo[i].w << std::endl;
			// viewport transform on each vertex
			ptsHomo[i].x = ptsHomo[i].x * (width/2.0f) + ptsHomo[i].w * (width/2.0f);
			ptsHomo[i].y = ptsHomo[i].y * (height/2.0f) + ptsHomo[i].w * (height/2.0f);
			// std::cout << ptsHomo[i].x << " " << ptsHomo[i].y << std::endl;
			minX = std::min(minX, ptsHomo[i].x);
			maxX = std::max(maxX, ptsHomo[i].x);
			minY = std::min(minY, ptsHomo[i].y);
			maxY = std::max(maxY, ptsHomo[i].y);
		}
		glm::vec4 segs[3] = {
			ptsHomo[1] - ptsHomo[0],
			ptsHomo[2] - ptsHomo[1],
			ptsHomo[0] - ptsHomo[2]
		};

		// rasterization
		for(int x = floor(minX); x < ceil(maxX); x += 1){
			for(int y = floor(minY); y < ceil(maxY); y += 1){
				if(x + 0.5f < 0 || x + 0.5f >= width){
					continue;
				}
				if(y + 0.5f < 0 || y + 0.5f >= height){
					continue;
				}
				float x_ctr = x + 0.5f;
				float y_ctr = y + 0.5f;
				int mn_sign = 1e9; int mx_sign = -1e9;
				for(int j = 0; j < 3; j++){
					float to_pt[2] = {x_ctr - ptsHomo[(j + 1) % 3].x, y_ctr - ptsHomo[(j + 1) % 3].y};
					float cross = segs[j].x * to_pt[1] - segs[j].y * to_pt[0];
					cross = glm::sign(cross);
					mn_sign = std::min(mn_sign, (int) cross);
					mx_sign = std::max(mx_sign, (int) cross);
				}
				if(abs(mx_sign - mn_sign) > 1){
					// not inside
					continue;
				}
				glm::vec3 curr_px(x_ctr, y_ctr, 0);
				glm::vec3 bary = getBarycentricConstants(curr_px, ptsHomo);
				float z_avg = bary.x * ptsHomo[0].z + bary.y * ptsHomo[1].z + bary.z * ptsHomo[2].z;
				// z_avg /= 3.0f;
				if(at(depth, y, x) > z_avg){
					// closer to camera
					set(depth, y, x, -1, z_avg);
					// show this pixel
					glm::vec3 px_col;
					if(tex_mode == -1){
						// in colored mode
						px_col = c[0] * bary.x + c[1] * bary.y + c[2] * bary.z;
					} else {
						// in textured mode
						float tex_x, tex_y;
						getTexCoords(curr_px, ptsHomo, z_inv, Q_sca, tex_x, tex_y);
						size_t img_x = floor(tex_x * tex_width[0]);
						size_t img_y = floor(tex_y * tex_height[0]);
						size_t img_x_close = round(tex_x * tex_width[0]);
						size_t img_y_close = round(tex_y * tex_height[0]);
						// printf("%d %d %.5f %.5f\n", x, y, tex_x, tex_y);
						if(tex_mode == 0){
							// nn
							px_col.x = at_tex(textures[0], img_y, img_x, 0, 0);
							px_col.y = at_tex(textures[0], img_y, img_x, 0, 1);
							px_col.z = at_tex(textures[0], img_y, img_x, 0, 2);
						} else if (tex_mode == 1){
							// bilinear
							px_col = bilinear_interp(0, img_x_close, img_y_close, tex_x, tex_y);
						} else {
							// mipmap
							// get tex coords at (x, y), (x + 1, y), (x, y + 1)
							float tex_x_curr, tex_y_curr;
							float tex_x_up, tex_y_up;
							float tex_x_right, tex_y_right;
							glm::vec3 px_right(curr_px.x + 1, curr_px.y, curr_px.z);
							glm::vec3 px_up(curr_px.x, curr_px.y + 1, curr_px.z);
							// printf("passed 1\n");
							getTexCoords(curr_px, ptsHomo, z_inv, Q_sca, tex_x_curr, tex_y_curr);
							getTexCoords(px_right, ptsHomo, z_inv, Q_sca, tex_x_right, tex_y_right);
							getTexCoords(px_up, ptsHomo, z_inv, Q_sca, tex_x_up, tex_y_up);
							// printf("passed 2\n");
							float first = sqrt(pow((tex_x_right - tex_x_curr) * tex_width[0], 2) + pow((tex_y_right - tex_y_curr) * tex_height[0], 2)); // du/dx + dv/dx
							float second = sqrt(pow((tex_x_up - tex_x_curr) * tex_width[0], 2) + pow((tex_y_up - tex_y_curr) * tex_height[0], 2)); // du/dy + dv/dy
							float L = max(first, second);
							float D = log2(L);
							int D_lo = floor(D);
							int D_hi = ceil(D);
							D_lo = CLAMP(D_lo, 0, textures.size() - 1);
							D_hi = CLAMP(D_hi, 0, textures.size() - 1);
							img_x_close = round(tex_x * tex_width[0]);
							img_y_close = round(tex_y * tex_height[0]);
							// printf("%d %d\n", tex_width[0], tex_width[D_lo]);
							glm::vec3 col_first = bilinear_interp(D_lo, img_x_close, img_y_close, tex_x, tex_y);
							img_x_close = round(tex_x * tex_width[0]);
							img_y_close = round(tex_y * tex_height[0]);
							glm::vec3 col_second = bilinear_interp(D_hi, img_x_close, img_y_close, tex_x, tex_y);
							float w_1 = 1.0f - (D - D_lo);
							float w_2 = 1.0f - w_1;
							px_col = w_1 * col_first + w_2 * col_second;
						}
					}
					set(color, y, x, 0, px_col.x);
					set(color, y, x, 1, px_col.y);
					set(color, y, x, 2, px_col.z);
				}
			}
		}
}