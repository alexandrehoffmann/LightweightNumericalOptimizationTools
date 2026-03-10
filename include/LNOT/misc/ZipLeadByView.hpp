#ifndef LNOT_MISC_ZIP_LEAD_VIEW_BY_HPP
#define LNOT_MISC_ZIP_LEAD_VIEW_BY_HPP

#include <ranges>

namespace LNOT::misc
{
	
template<class... Views>
concept CZipableViews = (std::ranges::input_range<Views> and ...) and (std::ranges::view<Views> and ...) and requires
{
	typename std::common_type_t<std::ranges::range_difference_t<Views>...>;
};

template<class... Ranges>
concept CZipableRanges = CZipableViews<std::views::all_t<Ranges>...>;

template<std::size_t I, class... Views>
class ZipLeadByView : public std::ranges::view_interface<ZipLeadByView<I, Views...>>
{
	static_assert(CZipableViews<Views...>);
	static_assert(I < sizeof...(Views));
public:
	using LeadingView = std::tuple_element_t<I, std::tuple<Views...> >;
    
	static constexpr bool hasConst = (std::ranges::range<Views> and ...);
    
	template<bool isConst>
	using Sentinel = std::ranges::sentinel_t< std::conditional_t<isConst, const LeadingView, LeadingView> >;

	template<bool isConst>
	class Iterator
	{
		template<std::ranges::view View>
		using MaybeConst = std::conditional_t<isConst, std::add_const_t<View>, View>;
		
		template<std::ranges::view View>
		using IteratorFor = std::ranges::iterator_t< MaybeConst<View> >;

		template<std::ranges::view View>
		using SentinelFor = std::ranges::sentinel_t< MaybeConst<View> >;

		//TODO: switch to iter_const_reference_t whenever moving to C++23
		template<std::ranges::view View>
		using ReferenceFor = std::common_reference_t<
			const std::iter_value_t<IteratorFor<View>>&&,
			std::iter_reference_t<IteratorFor<View>>>;

		template<std::ranges::view View>
		using ValueFor = std::iter_value_t< IteratorFor<View> >;

		template<std::ranges::view View>
		using DifferenceFor = std::iter_difference_t< IteratorFor<View> >;
		
	public:
		using value_type      = std::tuple<ValueFor<Views>...>;
		using reference       = std::tuple<ReferenceFor<Views>...>;
		using difference_type = std::common_type_t<DifferenceFor<Views>...>;

		static constexpr bool isInput         = (std::input_iterator< IteratorFor<Views> >         and ...);
		static constexpr bool isForward       = (std::forward_iterator< IteratorFor<Views> >       and ...);
		static constexpr bool isBidirectional = (std::bidirectional_iterator< IteratorFor<Views> > and ...);
		static constexpr bool isRandomAccess  = (std::random_access_iterator< IteratorFor<Views> > and ...);
		static constexpr bool isContiguous    = (std::contiguous_iterator< IteratorFor<Views> >    and ...);

		static_assert(isInput);

		using iterator_category = std::conditional_t<isContiguous, 
		    std::contiguous_iterator_tag, 
		    std::conditional_t<isRandomAccess, 
				std::random_access_iterator_tag,
				std::conditional_t<isBidirectional,
					std::bidirectional_iterator_tag,
					std::conditional_t<isForward,
						std::forward_iterator_tag,
						std::input_iterator_tag>
					>
				>
		>;
		    
		constexpr explicit Iterator(const IteratorFor<Views>... iterators) : m_iterators(iterators...) {}
		
		//=== input iterator
		
		constexpr reference operator*() const requires(isInput);
		
		constexpr Iterator& operator++() requires(isInput);
		
		constexpr Iterator operator++(int) requires(isInput) { Iterator ret(*this); ++(*this); return ret; }
		
		friend constexpr bool operator==(const Iterator& it, const SentinelFor<LeadingView> bound) requires(isInput) { return std::get<I>(it.m_iterators) == bound; }
		
		//==== forward iterator

		friend constexpr bool operator==(const Iterator& lhs, const Iterator& rhs) requires(isForward) { return lhs.m_iterators == rhs.m_iterators; }
		
		//==== bidirectional iterator

		constexpr Iterator& operator--() requires(isBidirectional);
		
		constexpr Iterator operator--(int) requires(isBidirectional) { Iterator ret(*this); --(*this); return ret; }
		
		//==== random access iterator
		
		constexpr Iterator& operator+=(const difference_type n) requires(isRandomAccess);
		
		constexpr Iterator& operator-=(const difference_type n) requires(isRandomAccess);
		
		constexpr reference operator[](const difference_type n) const requires(isRandomAccess);
		
		friend constexpr Iterator operator+(const Iterator& it, const difference_type n) requires(isRandomAccess) { Iterator ret(it); return it += n; }
		
		friend constexpr Iterator operator+(const difference_type n, const Iterator& it) requires(isRandomAccess) { return it + n; }
		
		friend constexpr Iterator operator-(const Iterator& it, const difference_type n) requires(isRandomAccess) { Iterator ret(it); return it -= n; }
		
		friend constexpr difference_type operator-(const Iterator& lhs, const Iterator& rhs) requires(isRandomAccess) { return std::get<I>(lhs.m_iterators) - std::get<I>(rhs.m_iterators); }
		
		friend constexpr auto operator<=>(const Iterator& lhs, const Iterator& rhs) requires(isRandomAccess) { return lhs.m_iterators <=> rhs.m_iterators; } 

	private:
		std::tuple<IteratorFor<Views>...> m_iterators;
	};
	
	constexpr explicit ZipLeadByView(Views&&... views) : m_views(std::forward<Views>(views)...) {} 
	
	constexpr Iterator<false> begin();
	
	constexpr Iterator<true> begin() requires(hasConst);
	
	constexpr Sentinel<false> end() { return std::ranges::end(std::get<I>(m_views)); } 
	
	constexpr Sentinel<true> end() requires(hasConst) { return std::ranges::end(std::get<I>(m_views)); }
	
	constexpr std::unsigned_integral auto size() const requires(std::ranges::sized_range<LeadingView>) { return std::get<I>(m_views).size(); }
private:
	std::tuple<Views...> m_views;
};

template<std::size_t I, std::ranges::range... Ranges>
constexpr std::ranges::view auto zipLeadBy(Ranges&&... ranges) requires(CZipableRanges<Ranges...>) { return ZipLeadByView<I, std::views::all_t<Ranges>...>(std::views::all(ranges)...); }

} // namespace LNOt::misc

#include <LNOT/misc/ZipLeadByView_impl.hpp>

#endif // LNOT_MISC_ZIP_LEAD_VIEW_BY_HPP
