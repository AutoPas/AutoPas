/**
 * @file StorageClass.h
 * @author J. Körner
 * @date 19.04.2021
 */
#pragma once

template<class domainId>
class StorageClass {
	public:
		typedef DomainId domainId;

		virtual void getNeighboursOfSubdomain(DomainId subdomainIndex, std::vector<DomainId> neighbours) = 0;

	protected
		StorageClass() = default;
		~StorageClass() = default;
}


