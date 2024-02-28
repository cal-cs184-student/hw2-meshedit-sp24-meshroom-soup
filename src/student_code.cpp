#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    vector<Vector2D> cur_level = points;
    vector<Vector2D> next_level;
    for (int i = 0; i < cur_level.size() - 1; i++) {
      next_level.push_back((1 - t) * cur_level[i] + t * cur_level[i+1]);
    }
    return next_level;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    vector<Vector3D> cur_level = points;
    vector<Vector3D> next_level;
    for (int i = 0; i < cur_level.size() - 1; i++) {
      next_level.push_back((1 - t) * cur_level[i] + t * cur_level[i+1]);
    }
    return next_level;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    vector<Vector3D> cur_level = points;
    while (cur_level.size() != 1) {
      cur_level = evaluateStep(cur_level, t);
    }
    return cur_level[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
    vector<Vector3D> col_points;
    for (int i = 0; i < controlPoints.size(); i++) {
      col_points.push_back(evaluate1D(controlPoints[i], u));
    }
    return evaluate1D(col_points, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
    HalfedgeCIter h1 = halfedge();
    HalfedgeCIter h2 = h1->next();
    HalfedgeCIter h3 = h2->next();
    Vector3D weighted_norm_sum = Vector3D();
    do {
      // get area-weighted norm and add to total
      Vector3D AB = h2->vertex()->position - position;
      Vector3D BC = h3->vertex()->position - h2->vertex()->position;
      weighted_norm_sum += cross(AB, BC) / 2;

      // move h1, h2, h3 for next adjacent triangle
      h1 = h3->twin();
      h2 = h1->next();
      h3 = h2->next();
    } while (h1 != halfedge());
    return weighted_norm_sum / weighted_norm_sum.norm();
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (e0->isBoundary()) return e0;

    HalfedgeIter h1 = e0->halfedge();
    HalfedgeIter h2 = h1->twin();
    HalfedgeIter h3 = h1->next();
    HalfedgeIter h4 = h2->next();
    HalfedgeIter h5 = h3->next();
    HalfedgeIter h6 = h4->next();
    VertexIter v1 = h1->vertex();
    VertexIter v2 = h2->vertex();
    FaceIter f1 = h1->face();
    FaceIter f2 = h2->face();

    h1->setNeighbors(h6, h2, h5->vertex(), e0, f1);
    h2->setNeighbors(h5, h1, h6->vertex(), e0, f2);
    h3->setNeighbors(h1, h3->twin(), h3->vertex(), h3->edge(), f1);
    h4->setNeighbors(h2, h4->twin(), h4->vertex(), h4->edge(), f2);
    h5->setNeighbors(h4, h5->twin(), h5->vertex(), h5->edge(), f2);
    h6->setNeighbors(h3, h6->twin(), h6->vertex(), h6->edge(), f1);
    v1->halfedge() = h4;
    v2->halfedge() = h3;
    f1->halfedge() = h1;
    f2->halfedge() = h2;

    return EdgeIter();
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
    if (e0->isBoundary()) return VertexIter();

    HalfedgeIter h1 = e0->halfedge();
    HalfedgeIter h2 = h1->twin();
    HalfedgeIter h3 = h1->next();
    HalfedgeIter h4 = h2->next();
    HalfedgeIter h5 = h3->next();
    HalfedgeIter h6 = h4->next();
    VertexIter v1 = h1->vertex();
    VertexIter v2 = h2->vertex();
    FaceIter f1 = h1->face();
    FaceIter f2 = h2->face();

    VertexIter m = newVertex();
    m->position = (v1->position + v2->position) / 2;
    m->halfedge() = h1;
    m->isNew = true;
    // m->newPosition = m->position;
    // m->centroid = v1->centroid;
    // m->quadric = v1->quadric;

    HalfedgeIter h7 = newHalfedge();
    HalfedgeIter h8 = newHalfedge();
    HalfedgeIter h9 = newHalfedge();
    HalfedgeIter h10 = newHalfedge();
    HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    EdgeIter e1 = newEdge();
    EdgeIter e2 = newEdge();
    EdgeIter e3 = newEdge();
    FaceIter f3 = newFace();
    FaceIter f4 = newFace();

    // updating halfedges
    h1->setNeighbors(h1->next(), h2, m, e0, f1);
    h2->setNeighbors(h7, h1, h2->vertex(), e0, f2);
    h3->setNeighbors(h12, h3->twin(), h3->vertex(), h3->edge(), f1);
    h4->setNeighbors(h8, h4->twin(), h4->vertex(), h4->edge(), f4);
    h5->setNeighbors(h10, h5->twin(), h5->vertex(), h5->edge(), f3);
    h6->setNeighbors(h2, h6->twin(), h6->vertex(), h6->edge(), f2);
    h7->setNeighbors(h6, h8, m, e1, f2);
    h8->setNeighbors(h9, h7, h6->vertex(), e1, f4);
    h9->setNeighbors(h4, h10, m, e2, f4);
    h10->setNeighbors(h11, h9, h4->vertex(), e2, f3);
    h11->setNeighbors(h5, h12, m, e3, f3);
    h12->setNeighbors(h1, h11, h5->vertex(), e3, f1);

    // updating vertices
    h3->vertex()->halfedge() = h3;
    h5->vertex()->halfedge() = h5;
    h4->vertex()->halfedge() = h4;
    h6->vertex()->halfedge() = h6;

    // updating edges
    e0->halfedge() = h1;
    e1->halfedge() = h7;
    e2->halfedge() = h9;
    e3->halfedge() = h11;
    e0->isNew = false;
    e1->isNew = true;
    e2->isNew = true;
    e3->isNew = true;
    e0->isBlue = false;
    e1->isBlue = true;
    e2->isBlue = false;
    e3->isBlue = true;

    // updating faces
    f1->halfedge() = h1;
    f2->halfedge() = h2;
    f3->halfedge() = h5;
    f4->halfedge() = h4;

    return m;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      // iterate over all neighboring vertices to get new positions of old vertices
      HalfedgeCIter h = v->halfedge(); 
      int n = v->degree();
      float u = (n == 3) ? 3.0 / 16 : 3.0 / (8 * n);
      Vector3D orig_nei_pos_sum = Vector3D();
      do {
        HalfedgeCIter h_twin = h->twin();
        VertexCIter v_nei = h_twin->vertex();
        orig_nei_pos_sum += v_nei->position;
        h = h_twin->next();
      } while(h != v->halfedge()); 
      v->newPosition = (1 - n * u) * v->position + u * orig_nei_pos_sum;
      v->isNew = false;
    }
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      Vector3D v1_pos = e->halfedge()->vertex()->position;
      Vector3D v2_pos = e->halfedge()->twin()->vertex()->position;
      Vector3D v3_pos = e->halfedge()->next()->next()->vertex()->position;
      Vector3D v4_pos = e->halfedge()->twin()->next()->next()->vertex()->position;
      e->newPosition = 3.0 / 8 * (v1_pos + v2_pos) + 1.0 / 8 * (v3_pos + v4_pos);
      e->isNew = false;
      e->isBlue = false;
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      if (!e->isNew && !e->isBlue) {
        VertexIter v = mesh.splitEdge(e);
        v->newPosition = e->newPosition;
      }
    }
    
    // 4. Flip any new edge that connects an old and new vertex.
    for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
      VertexIter v1 = e->halfedge()->vertex();
      VertexIter v2 = e->halfedge()->twin()->vertex();
      if (v1->isNew != v2->isNew && e->isBlue) {
        mesh.flipEdge(e);
      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
      v->position = v->newPosition;
    }
  }
}
