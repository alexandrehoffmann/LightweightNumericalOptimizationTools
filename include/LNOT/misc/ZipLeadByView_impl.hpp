#ifndef LNOT_MISC_ZIP_LEAD_VIEW_BY_IMPL_HPP
#define LNOT_MISC_ZIP_LEAD_VIEW_BY_IMPL_HPP

namespace LNOT::misc
{

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator*() const -> reference requires(isInput)
{
	return std::apply([](const IteratorFor<Views>&... iterators) -> reference
	{
		return reference((*iterators)...);
	}, m_iterators);
}

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator++() -> Iterator& requires(isInput)
{
	std::apply([](IteratorFor<Views>&... iterators) -> void
	{
		(++iterators, ...);
	}, m_iterators);

	return *this;
}

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator--() -> Iterator& requires(isBidirectional) 
{
	std::apply([](IteratorFor<Views>&... iterators) -> void
	{
		(--iterators, ...);
	}, m_iterators);

	return *this;
}

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator+=(const difference_type n) -> Iterator& requires(isRandomAccess)
{
	std::apply([n](IteratorFor<Views>&... iterators) -> void
	{
		((iterators += n), ...);
	}, m_iterators);

	return *this;
}

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator-=(const difference_type n) -> Iterator& requires(isRandomAccess)
{
	std::apply([n](IteratorFor<Views>&... iterators) -> void
	{
		((iterators -= n), ...);
	}, m_iterators);

	return *this;
}

template<std::size_t I, class... Views> template<bool isConst>
constexpr auto ZipLeadByView<I, Views...>::Iterator<isConst>::operator[](const difference_type n) const -> reference requires(isRandomAccess)
{
	return std::apply([n](const IteratorFor<Views>&... iterators) -> reference
	{
		return reference((iterators[n])...);
	}, m_iterators);
}

template<std::size_t I, class... Views> 
constexpr auto ZipLeadByView<I, Views...>::begin() -> Iterator<false>
{
	return std::apply([](Views&... views) -> Iterator<false>
	{
		return Iterator<false>(std::ranges::begin(views)...);
	}, m_views);
}

template<std::size_t I, class... Views> 
constexpr auto ZipLeadByView<I, Views...>::begin() -> Iterator<true> requires(hasConst)
{
	return std::apply([](Views&... views) -> Iterator<true>
	{
		return Iterator<true>(std::ranges::begin(views)...);
	}, m_views);
}

} // namespace LNOt::misc

#endif // LNOT_MISC_ZIP_LEAD_VIEW_BY_IMPL_HPP
