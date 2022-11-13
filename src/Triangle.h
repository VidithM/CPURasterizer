#pragma once

#include <stdlib.h>
#include <math.h>
#include <vector>

#include <glm/glm.hpp>


class Triangle {
	private:
		glm::vec3 v[3];		// Triangle vertices
		glm::vec3 c[3];		// Vertex color
		glm::vec2 t[3];		// Texture coordinates

	public:

		// Default constructor
		Triangle();

		// Constructor without texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2);

		// Constructor with texture coordinates
		Triangle(glm::vec3& v0, glm::vec3& v1, glm::vec3& v2, glm::vec2& t0, glm::vec2& t1, glm::vec2& t2);

		// Rendering the triangle using OpenGL
		void RenderOpenGL(glm::mat4 &modelViewMatrix, glm::mat4 &projectionMatrix, bool textureMode);

		void setColor(int vert, glm::vec3& color);

		glm::vec3 getVert(int vert);

		// glm::vec3 getBarycentricConstants(glm::vec3& pos);

		// Rendering the triangle using CPU
		void RenderCPU(
			glm::mat4 &modelViewMatrix, 
			glm::mat4 &projectionMatrix,
			float height,
			float width,
			float *color,
			std::vector<float*> textures,
			std::vector<size_t> tex_height,
			std::vector<size_t> tex_width,
			int tex_mode,
			float *depth
		);
};
