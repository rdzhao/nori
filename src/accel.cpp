/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/accel.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    /* Nothing to do here for now */
	std::list<int> tris;
	for (int i = 0; i < m_mesh->getTriangleCount(); ++i)
		tris.push_back(i);
	cout << "Total Triangle Num: " << m_mesh->getTriangleCount() << endl;
	root = buildOctNode(m_mesh->getBoundingBox(), tris, 0);

	cout << "Build Complete .. "<< endl;
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    ///* Brute force search through all triangles */
    //for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //    float u, v, t;
    //    if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //        /* An intersection was found! Can terminate
    //           immediately if this is a shadow ray query */
    //        if (shadowRay)
    //            return true;
    //        ray.maxt = its.t = t;
    //        its.uv = Point2f(u, v);
    //        its.mesh = m_mesh;
    //        f = idx;
    //        foundIntersection = true;
    //    }
    //}

	/*Oct Tree recursive search*/
	//cout << "Before intersection test ..."<< endl;
	rayIntersectRecursive(root, ray, its, f, foundIntersection, shadowRay);
	//cout << "After test " << foundIntersection << " " << f << endl;

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }
	//cout << "Intersection Process done! "<< endl;
    return foundIntersection;
}

OctNodeIter Accel::buildOctNode(BoundingBox3f bb, std::list<int> tris, int depth)
{
	if (tris.empty())
		return nodes.end();

	if (tris.size() < 10 || depth > MAXDEPTH){
		OctNode octNode;
		nodes.push_back(octNode);
		OctNodeIter otIter = --nodes.end();
		
		otIter->bbox = bb;
		otIter->triangles = tris;

		return otIter;
	} 
	else{
		std::vector<BoundingBox3f> subBoxes;
		std::vector<std::list<int>> subLists;
		subBoxes.resize(8);
		subLists.resize(8);

		for (int j = 0; j < 8; ++j) {
			Vector3f center, corner;
			center = bb.getCenter();
			corner = bb.getCorner(j);

			BoundingBox3f subBox(center);
			subBox.expandBy(corner);
			subBoxes[j] = subBox;
		}

		for (auto iter = tris.begin(); iter != tris.end(); ++iter) {
			for (int j = 0; j < 8; ++j){
				if (subBoxes[j].overlaps(m_mesh->getBoundingBox(*iter))){
					subLists[j].push_back(*iter);
				}
			}
		}

		OctNode octNode;
		nodes.push_back(octNode);
		OctNodeIter otIter = --nodes.end();

		otIter->bbox = bb;
		for (int i = 0; i < 8; ++i) {
			otIter->children.push_back(buildOctNode(subBoxes[i], subLists[i], depth + 1));
		}

		subBoxes.clear();
		subLists.clear();

		return otIter;
	}
}

void Accel::rayIntersectRecursive(OctNodeIter node, Ray3f &ray, Intersection &its, uint32_t& f, bool& found, bool shadowRay) const
{
	if (node->triangles.empty()){ 
		// reorder here
		//std::sort(node->children.begin(), node->children.end(), Compare(&ray));

		for (int i = 0; i < 8; ++i){
			if (node->children[i] != nodes.end() && node->children[i]->bbox.rayIntersect(ray))
				rayIntersectRecursive(node->children[i], ray, its, f, found, shadowRay);
			//if (found)
				//break;
		}
	} // not leaf
	else {
		for (auto iter = node->triangles.begin(); iter != node->triangles.end(); ++iter){
			float u, v, t;
			if (m_mesh->rayIntersect(*iter, ray, u, v, t)) {
				/* An intersection was found! Can terminate
				   immediately if this is a shadow ray query */
				if (shadowRay)
					found = true;
				ray.maxt = its.t = t;
				its.uv = Point2f(u, v);
				its.mesh = m_mesh;
				f = *iter;
				found = true;
			}
		}
	} // leaf
}

NORI_NAMESPACE_END

