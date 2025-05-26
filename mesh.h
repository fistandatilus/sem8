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
    unsigned int get_matrix_from_mesh(const unsigned int i, const unsigned int j) const
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
                return 3 * width + (j - width) * width + (i - 2 * width);
            if (j == 2 * width && i < 3 * width)
                return 4 * width * width + i;
            return 4 * width * width + 3 * width + j;
        }

        return 0;
    }
    void get_mesh_from_matrix(unsigned int i_matrix, unsigned int &i, unsigned int &j) const
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

    mesh(const unsigned int width) : width(width)
    {};

private:
    unsigned int width{}; //amount of subdivisions of [0,1]
};

#endif //MESH_H
