#ifndef LNOT_CIRCULAR_BUFFER_IMPL_HPP
#define LNOT_CIRCULAR_BUFFER_IMPL_HPP

#include <LNOT/Containers/CircularBuffer.hpp>

namespace LNOT
{
	
//// explicit template instanciations ////

extern template class CircularBuffer<float>;
extern template class CircularBuffer<double>;

//// method implementations ////
	
template<typename T, typename Allocator>
CircularBuffer<T, Allocator>::CircularBuffer(const size_type size, const allocator& alloc) : 
	m_allocator  (alloc),
	m_data_begin (m_allocator.allocate (size)),
	m_data_end   (m_data_begin + size), m_range_end (m_data_begin),
	m_allFilled  (false)
{
}

template<typename T, typename Allocator>
CircularBuffer<T, Allocator>::CircularBuffer(const CircularBuffer &other) : 
	m_data_begin (m_allocator.allocate (other.size ())),
	m_data_end   (m_data_begin + other.size ()),
	m_range_end  (m_data_begin + (other.m_range_end - other.m_data_begin)),
	m_allFilled  (other.m_allFilled)
{
	std::copy(other.m_data_begin, other.m_data_end, m_data_begin);
}

template<typename T, typename Allocator>
CircularBuffer<T, Allocator>::~CircularBuffer() 
{ 
	for (pointer it=m_data_begin; it!=m_data_end; ++it) {  std::destroy_at(it); } 
	m_allocator.deallocate(m_data_begin, size()); 
}

template<typename T, typename Allocator> template<typename... Args>
auto CircularBuffer<T, Allocator>::push(Args&&...args) -> size_type
{
	if (m_range_end == m_data_end)
	{
		m_range_end = m_data_begin;
		m_allFilled = true;
	}
	std::construct_at(m_range_end, std::forward<Args>(args)...);
	++m_range_end;
	
	return size_type(m_range_end - 1 - m_data_begin);
}

template<typename T, typename Allocator> template<typename Function> 
Function CircularBuffer<T, Allocator>::foreach(Function f) requires (IsTraverserFunction<Function>::value)
{
	if (m_allFilled and m_range_end != m_data_end)
	{
		bool isDone = false;
		for (pointer it=m_range_end; not isDone; ++it)
		{
			if (it == m_data_end) { it = m_data_begin; }
			f(size_type(it - m_data_begin), *it);
			isDone = (std::next(it) == m_range_end);
		}
	}
	else
	{
		for (pointer it=m_data_begin; it!=m_range_end; ++it)
		{
			f(size_type(it - m_data_begin), *it);
		}
	}
	return f;
}

template<typename T, typename Allocator> template<typename Function> 
Function CircularBuffer<T, Allocator>::foreach(Function f) const requires (IsConstTraverserFunction<Function>::value)
{
	if (m_allFilled and m_range_end != m_data_end)
	{
		bool isDone = false;
		for (pointer it=m_range_end; not isDone; ++it)
		{
			if (it == m_data_end) { it = m_data_begin; }
			f(size_type(it - m_data_begin), *it);
			isDone = (std::next(it) == m_range_end);
		}
	}
	else
	{
		for (pointer it=m_data_begin; it!=m_range_end; ++it)
		{
			f(size_type(it - m_data_begin), *it);
		}
	}
	return f;
}

template<typename T, typename Allocator> template<typename Function> 
Function CircularBuffer<T, Allocator>::reverseForeach(Function f) requires (IsTraverserFunction<Function>::value)
{  
    const pointer data_rbegin = std::prev(m_data_end);
    const pointer data_rend = std::prev(m_data_begin);
    if (m_allFilled and m_range_end != m_data_end)
    {
        bool isDone = false;
        for (pointer it=std::prev(m_range_end); not isDone; --it)
        {
            if (it == data_rend) { it = data_rbegin; }
						f(size_type(it - m_data_begin), *it);
            isDone = (it == m_range_end);
        }
    }
    else
    {
        for (pointer it=std::prev(m_range_end); it!=data_rend; --it)
        {
					f(size_type(it - m_data_begin), *it);
        }
    }
    return f;
}

template<typename T, typename Allocator> template<typename Function> 
Function CircularBuffer<T, Allocator>::reverseForeach(Function f) const requires (IsConstTraverserFunction<Function>::value)
{  
    const const_pointer data_rbegin = std::prev(m_data_end);
    const const_pointer data_rend = std::prev(m_data_begin);
    if (m_allFilled and m_range_end != m_data_end)
    {
        bool isDone = false;
        for (const_pointer it=std::prev(m_range_end); not isDone; --it)
        {
            if (it == data_rend) { it = data_rbegin; }
						f(it - m_data_begin, *it);
            isDone = (it == m_range_end);
        }
    }
    else
    {
        for (const_pointer it=std::prev(m_range_end); it!=data_rend; --it)
        {
					f(it - m_data_begin, *it);
        }
    }
    return f;
}


} // namespace LNOT

#endif // LNOT_CIRCULAR_BUFFER_IMPL_HPP
