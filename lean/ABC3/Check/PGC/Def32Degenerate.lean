import ABC3.Skeleton.PGC.Section3

/-!
# [pGC] Definition 3.2 の**旧**形は空虚だった——`I = ∅` が許されていた

原文 (pGC p.6):

> We shall call the E[Γ_K]-module V uniformizing if the restriction of ρ_V to some open
> subgroup I of U_K (⊆ Γ^a_K^b) is the morphism I → E× induced by restricting some morphism
> of fields K → E to I ⊆ U_K ⊆ K.

★注目すべきは "open **subgroup**" ——「開集合」ではなく「開**部分群**」である。

`Skeleton/PGC/Section3.lean::IsUniformizing` の旧形は `I : Set K.carrier` に
`IsOpen I` と `I ⊆ U_K` しか課しておらず、**「部分群」の条件が落ちていた**。
すると `I = ∅` が許され、`∀ x ∈ I, …` は空虚に真になる。`E := K.carrier`・
`ι := RingHom.id` と取れば

**`IsUniformizing` はどんな `ρ` でも常に成り立つ**

——定義が内容を失う(`isUniformizingOld_trivial`、`sorry` 無し)。
そして Corollary 3.3 の `↔` も自明に真になってしまう。

★これは `Check/PGC/InertiaDegeneracy.lean`(`I_K` を自由データに置くと退化)・
`Check/PGC/Theorem42Degenerate.lean`(`Φ` が自由な関数だと偽)に続く3例目——
**落とした条件は、主張を偽にするか自明にするかのどちらかになる**。

## 修理

`Skeleton/PGC/Section3.lean::IsUniformizing` に原文どおり
「`I` は `U_K` の部分群」(`1 ∈ I`・積で閉じる・逆元で閉じる)を課した。
`Lemma 4.1`(`Section4.lean`)が `1 ∈ I` を課していたのと同じ形である。

修理後は定義が**識別力を持つ**: 自明な表現 `ρ = 1` は uniformizing に
**ならない**(`not_isUniformizing_one`)——`I` が開で `1 ∈ I` なら
`1` の近傍を含み、体準同型 `ι` は単射なので `ι x = 1` は `x = 1` を強いる。
`K` の位相は非離散なのでそれは不可能。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-- **旧** Definition 3.2——「部分群」の条件が無い形(記録用に写したもの)。 -/
def IsUniformizingOld (K : PAdicLocalField p) (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) (ρ : K.absGal →* Eˣ) : Prop :=
  ∃ (I : Set K.carrier) (hIU : I ⊆ {x : K.carrier | ‖x‖ = 1}) (_hopen : IsOpen I)
    (ι : K.carrier →+* E), ∀ x (hx : x ∈ I), ((ρ (toGal ⟨x, hIU hx⟩) : Eˣ) : E) = ι x

/-- **★★★★★旧形は空虚に真**——`I = ∅`・`E := K`・`ι := id` で、
**どんな `ρ` でも** uniformizing になってしまう。 -/
theorem isUniformizingOld_trivial (K : PAdicLocalField p)
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) (ρ : K.absGal →* (K.carrier)ˣ) :
    IsUniformizingOld K K.carrier toGal ρ := by
  refine ⟨∅, Set.empty_subset _, isOpen_empty, RingHom.id K.carrier, ?_⟩
  intro x hx
  exact absurd hx (Set.notMem_empty x)

/-- **★★★★★修理後は識別力を持つ**——自明な表現 `ρ = 1` は uniformizing で**ない**。

`I` が開で `1 ∈ I` なら `1` のある近傍を含む。体準同型 `ι` は単射なので
`ι x = 1 = ι 1` は `x = 1` を強いるが、`K` の位相は非離散
(`NormedField.exists_norm_lt`)なのでそんな近傍は無い。 -/
theorem not_isUniformizing_one (K : PAdicLocalField p) (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) :
    ¬ IsUniformizing K E toGal 1 := by
  rintro ⟨I, hIU, hopen, hone, -, -, ι, hι⟩
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.mp hopen 1 hone
  obtain ⟨δ, hδpos, hδlt⟩ := NormedField.exists_norm_lt K.carrier hε
  have hmem : (1 + δ) ∈ I := by
    refine hball ?_
    simp only [Metric.mem_ball, dist_eq_norm]
    simpa using hδlt
  have h1 := hι (1 + δ) hmem
  have h2 := hι 1 hone
  simp only [MonoidHom.one_apply, Units.val_one] at h1 h2
  have h3 : ι (1 + δ) = ι 1 := by rw [← h1, ← h2]
  have h4 : (1 : K.carrier) + δ = 1 := ι.injective h3
  have h5 : δ = 0 := by linear_combination h4
  rw [h5, norm_zero] at hδpos
  exact lt_irrefl 0 hδpos

end ABC3.Check.PGC
