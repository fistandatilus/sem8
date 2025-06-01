#ifndef MESH_H
#define MESH_H

enum class point_type
{
    i,
    l,
    r,
    t,
    b,
    lt,
    rt,
    lb,
    rb,
    ilt,
    irb,
};

class mesh
{
public:
    unsigned int get_matrix_by_mesh(const unsigned int i, const unsigned int j) const
    {
        if (i < width)
        {
            if (j < width)
                return j * width + i;
            return 4 * width * width + i;
        }
        if (i < 2 * width)
        {
            if (j < 2 * width)
                return width * width + j * width + (i - width);
            return 4 * width * width + i;
        }
        if (i < 3 * width)
        {
            if (j >= width && j < 2 * width)
                return 3 * width * width + (j - width) * width + (i - 2 * width);
            if (j == 2 * width && i < 3 * width)
                return 4 * width * width + i;
        }
        return 4 * width * width + 3 * width + j;
    }
    void get_mesh_by_matrix(unsigned int i_matrix, unsigned int &i, unsigned int &j) const
    {
        if (i_matrix < width * width)
        {
            j = i_matrix / width;
            i = i_matrix % width;
            return;
        }
        i_matrix -= width * width;
        if (i_matrix < 2 * width * width)
        {
            j = i_matrix / width;
            i = width + i_matrix % width;
            return;
        }
        i_matrix -= 2 * width * width;
        if (i_matrix < width * width)
        {
            j = width + i_matrix / width;
            i = 2 * width + i_matrix % width;
            return;
        }
        i_matrix -= width * width;
        if (i_matrix < 3 * width)
        {
            i = i_matrix;
            j = 2 * width;
            if (i_matrix < width)
                j = width;
            return;
        }
        i_matrix -= 3 * width;
        i = 2 * width;
        if (i_matrix >= width)
            i = 3 * width;
        j = i_matrix;
    }

    point_type get_type_by_mesh(const unsigned int i, const unsigned int j) const
    {
        if (i == 0 && j == 0) 
            return point_type::lb;
        if (i == 3 * width && j == 2 * width)
            return point_type::rt;
        if ((i == 0 && j == width) ||
            (i == width && j == 2 * width))
            return point_type::lt;
        if ((i == 2 * width && j == 0) ||
            (i == 3 * width && j == width))
            return point_type::rb;
        if (i == width && j == width)
            return point_type::ilt;
        if (i == 2 * width && j == width)
            return point_type::irb;
        if ( i == 0 ||
            (i == width && j > width))
            return point_type::l;
        if ( i == 3 * width ||
            (i == 2 * width && j < width))
            return point_type::r;
        if ( j == 0 ||
            (i > 2 * width && j == width))
            return point_type::b;
        if ( j == 2 * width ||
            (i < width && j == width))
            return point_type::t;
        return point_type::i;
    }

    point_type get_type_by_matrix(const unsigned int i_matrix) const
    {
        unsigned int i;
        unsigned int j;
        get_mesh_by_matrix(i_matrix, i, j);

        return get_type_by_mesh(i, j);
    }

    mesh(const unsigned int width) : width(width), size(4 * width * width + 5 * width + 1)
    {};

    unsigned int width = 0; //amount of subdivisions of [0,1]
    unsigned int size = 0;
};

#endif //MESH_H
