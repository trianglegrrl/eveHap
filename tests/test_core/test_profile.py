"""Tests for evehap.core.profile module."""


from evehap.core.profile import AlleleObservation, AlleleProfile


class TestAlleleObservation:
    """Tests for AlleleObservation dataclass."""

    def test_create_observation(self) -> None:
        """Test creating an AlleleObservation."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 1.0},
            depth=100,
            quality=30.0,
            source="bam",
        )
        assert obs.position == 73
        assert obs.ref_allele == "A"
        assert obs.alleles == {"G": 1.0}
        assert obs.depth == 100
        assert obs.quality == 30.0
        assert obs.source == "bam"

    def test_major_allele(self) -> None:
        """Test major_allele property returns most frequent allele."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 0.95, "A": 0.05},
        )
        assert obs.major_allele == "G"

    def test_major_allele_tie_breaker(self) -> None:
        """Test major_allele with equal frequencies picks one consistently."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 0.5, "A": 0.5},
        )
        # Should return one of them consistently (alphabetically first)
        assert obs.major_allele in ["A", "G"]

    def test_is_variant_true(self) -> None:
        """Test is_variant returns True when major allele differs from ref."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 1.0},
        )
        assert obs.is_variant is True

    def test_is_variant_false(self) -> None:
        """Test is_variant returns False when major allele matches ref."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"A": 1.0},
        )
        assert obs.is_variant is False

    def test_is_heteroplasmic_true(self) -> None:
        """Test is_heteroplasmic returns True when multiple alleles above threshold."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 0.85, "A": 0.15},
        )
        assert obs.is_heteroplasmic is True

    def test_is_heteroplasmic_false(self) -> None:
        """Test is_heteroplasmic returns False when single allele dominant."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 0.98, "A": 0.02},
        )
        assert obs.is_heteroplasmic is False

    def test_is_heteroplasmic_single_allele(self) -> None:
        """Test is_heteroplasmic returns False for single allele."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 1.0},
        )
        assert obs.is_heteroplasmic is False

    def test_allele_frequency_present(self) -> None:
        """Test allele_frequency returns correct value for present allele."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 0.85, "A": 0.15},
        )
        assert obs.allele_frequency("G") == 0.85
        assert obs.allele_frequency("A") == 0.15

    def test_allele_frequency_absent(self) -> None:
        """Test allele_frequency returns 0.0 for absent allele."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 1.0},
        )
        assert obs.allele_frequency("T") == 0.0
        assert obs.allele_frequency("C") == 0.0

    def test_default_values(self) -> None:
        """Test default values for optional fields."""
        obs = AlleleObservation(
            position=73,
            ref_allele="A",
            alleles={"G": 1.0},
        )
        assert obs.depth is None
        assert obs.quality is None
        assert obs.source == "unknown"


class TestAlleleProfile:
    """Tests for AlleleProfile class."""

    def test_create_profile(self) -> None:
        """Test creating an AlleleProfile."""
        profile = AlleleProfile(sample_id="test_sample", source_format="bam")
        assert profile.sample_id == "test_sample"
        assert profile.source_format == "bam"
        assert len(profile.observations) == 0
        assert len(profile.covered_positions) == 0

    def test_add_observation(self) -> None:
        """Test adding an observation to a profile."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        obs = AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0})

        profile.add_observation(obs)

        assert 73 in profile.observations
        assert 73 in profile.covered_positions
        assert profile.observations[73] == obs

    def test_add_multiple_observations(self) -> None:
        """Test adding multiple observations."""
        profile = AlleleProfile(sample_id="test", source_format="bam")

        for pos in [73, 150, 263]:
            obs = AlleleObservation(position=pos, ref_allele="A", alleles={"G": 1.0})
            profile.add_observation(obs)

        assert len(profile.observations) == 3
        assert len(profile.covered_positions) == 3
        assert profile.covered_positions == {73, 150, 263}

    def test_get_allele_covered(self) -> None:
        """Test get_allele returns major allele for covered position."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        obs = AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0})
        profile.add_observation(obs)

        assert profile.get_allele(73) == "G"

    def test_get_allele_not_covered(self) -> None:
        """Test get_allele returns None for uncovered position."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert profile.get_allele(73) is None

    def test_get_observation_covered(self) -> None:
        """Test get_observation returns observation for covered position."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        obs = AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0})
        profile.add_observation(obs)

        assert profile.get_observation(73) == obs

    def test_get_observation_not_covered(self) -> None:
        """Test get_observation returns None for uncovered position."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert profile.get_observation(73) is None

    def test_is_derived_true(self) -> None:
        """Test is_derived returns True when position has derived allele."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        obs = AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0})
        profile.add_observation(obs)

        assert profile.is_derived(73, "G") is True

    def test_is_derived_false(self) -> None:
        """Test is_derived returns False when position doesn't have derived allele."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        obs = AlleleObservation(position=73, ref_allele="A", alleles={"A": 1.0})
        profile.add_observation(obs)

        assert profile.is_derived(73, "G") is False

    def test_is_derived_not_covered(self) -> None:
        """Test is_derived returns None when position not covered."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert profile.is_derived(73, "G") is None

    def test_coverage_fraction_empty(self) -> None:
        """Test coverage_fraction returns 0 for empty profile."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert profile.coverage_fraction() == 0.0

    def test_coverage_fraction_full(self) -> None:
        """Test coverage_fraction for full coverage."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        # mtDNA is 16569 bp
        for pos in range(1, 16570):
            obs = AlleleObservation(position=pos, ref_allele="A", alleles={"A": 1.0})
            profile.add_observation(obs)

        assert profile.coverage_fraction() == 1.0

    def test_coverage_fraction_partial(self) -> None:
        """Test coverage_fraction for partial coverage."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        # Cover ~50% of positions
        for pos in range(1, 8285):
            obs = AlleleObservation(position=pos, ref_allele="A", alleles={"A": 1.0})
            profile.add_observation(obs)

        # Should be approximately 50%
        frac = profile.coverage_fraction()
        assert 0.49 < frac < 0.51

    def test_mean_depth_no_data(self) -> None:
        """Test mean_depth returns None for empty profile."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        assert profile.mean_depth() is None

    def test_mean_depth_with_data(self) -> None:
        """Test mean_depth calculation."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.add_observation(
            AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0}, depth=100)
        )
        profile.add_observation(
            AlleleObservation(position=150, ref_allele="C", alleles={"T": 1.0}, depth=200)
        )
        profile.add_observation(
            AlleleObservation(position=263, ref_allele="A", alleles={"G": 1.0}, depth=300)
        )

        assert profile.mean_depth() == 200.0

    def test_mean_depth_with_missing_depth(self) -> None:
        """Test mean_depth skips observations without depth."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.add_observation(
            AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0}, depth=100)
        )
        profile.add_observation(
            AlleleObservation(position=150, ref_allele="C", alleles={"T": 1.0}, depth=None)
        )
        profile.add_observation(
            AlleleObservation(position=263, ref_allele="A", alleles={"G": 1.0}, depth=200)
        )

        assert profile.mean_depth() == 150.0

    def test_get_variants(self) -> None:
        """Test get_variants returns list of variants."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.add_observation(AlleleObservation(position=73, ref_allele="A", alleles={"G": 1.0}))
        profile.add_observation(
            AlleleObservation(position=150, ref_allele="C", alleles={"C": 1.0})  # No variant
        )
        profile.add_observation(AlleleObservation(position=263, ref_allele="A", alleles={"G": 1.0}))

        variants = profile.get_variants()
        assert len(variants) == 2
        assert (73, "A", "G") in variants
        assert (263, "A", "G") in variants

    def test_metadata(self) -> None:
        """Test metadata storage."""
        profile = AlleleProfile(sample_id="test", source_format="bam")
        profile.metadata["run_id"] = "ERR123456"
        profile.metadata["project"] = "AADR"

        assert profile.metadata["run_id"] == "ERR123456"
        assert profile.metadata["project"] == "AADR"
