/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GRID_DATAMANAGER_H_INCLUDED_
#define GRID_DATAMANAGER_H_INCLUDED_

#include <grid.h>

namespace grid
{
  struct octtree_piece_t
  {
    rect_t                         m_prct;
    rect_t                         m_ext_prct;
    int                            m_level;

    rect_t                         m_rect;
    rect_t                         m_ext_rect;
    rect_t                         m_domain_rect;

    octtree_piece_t(rect_t p,rect_t d,int l);

    std::string bn(const std::string& basename);

  };

  class data_manager_t
  {

  public:
    typedef boost::shared_ptr<octtree_piece_t> piece_ptr_t;
    typedef std::vector<piece_ptr_t>           piece_ptr_list_t;

  public:

    piece_ptr_list_t             m_pieces;

    cellid_t                     m_size;
    std::string                  m_filename;
    double                       m_simp_tresh;
    cell_fn_t                    m_f_range;
    cellid_t                     m_levels;

    std::string                  m_basename;
    int                          m_level_ct;

    data_manager_t
        ( std::string  filename,
          cellid_t     size,
          cellid_t     levels,
          double       simp_tresh
          );

    virtual ~data_manager_t ();

    void work();

    void createPieces();

    void split_dataset();

    void compute_subdomain_msgraphs ();

    void merge_up_subdomain_msgraphs ();

    void merge_down_subdomain_msgraphs ();

    void save_graphs ();

    void save_mfolds ();

    void destoryPieces();
  };
}

#endif
