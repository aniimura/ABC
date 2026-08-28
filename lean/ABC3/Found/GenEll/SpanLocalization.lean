/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HomogeneousCoordsArch
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★有限素点では座標の 1 つが他を割る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★これは何か —— `§9-940` の 1 本目の代数の部分

`§9-940` の**有限素点の整合**

    `∀ p Q, ∃ i y, (局所化した点が X_{s_i} を通る) ∧ (𝔞_Q = (x_i)) ∧ ((s_0/s_i)(y)·x_i = x_0)`

のうち、★★**後ろの 2 つは代数だけで出る**:

* `(𝓞_F)_Q` は**離散付値環**（`𝓞_F` は Dedekind、`Q ≠ ⊥`）
* 付値環では**割り切りが全順序**（`ValuationRing.dvd_total`）
* だから有限個の座標のうち**どれか 1 つが他を全部割る**——それを `x_j` とすれば
  `𝔞_Q = (x_j)` かつ `x_0/x_j ∈ (𝓞_F)_Q`

★★★これが「`v_Q(x_i)` が最小の `i` を取る」という古典的な操作の、
**付値を経由しない形**である。

## ★残っている段（明示）

★★残るのは**幾何の側の紐**である——ここで得た `j` が
`§9-939` の局所チャートの添字と**一致する**こと、すなわち

    `Spec (𝓞_F)_Q ⟶ X` が `X_{s_j}` を通り、そこでの `(s_0/s_j)` の値が `x_0/x_j` であること

★機構は `§9-913`（`ψ⁻¹(D₊(x_j)) = X_{s_j}`）と、座標が `ψ` から読まれていること
（`§9-941`）の合わせ技である。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★付値環では割り切りが全順序 -/

/-- ★★★★★★★★**付値環では有限族のどれか 1 つが他を全部割る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`ValuationRing.dvd_total`（割り切りが全順序）を有限集合の上で回すだけである。
★★これが「`v_Q(x_i)` が最小の `i` を取る」の**付値を経由しない形**である。 -/
theorem exists_dvd_all_of_valuationRing {S : Type} [CommRing S] [IsDomain S] [ValuationRing S]
    {ι : Type} [Fintype ι] [Nonempty ι] (y : ι → S) : ∃ j, ∀ k, y j ∣ y k := by
  classical
  have key : ∀ (s : Finset ι), s.Nonempty → ∃ j ∈ s, ∀ k ∈ s, y j ∣ y k := by
    intro s
    induction s using Finset.induction_on with
    | empty => intro h; exact absurd h (by simp)
    | insert a s ha ih =>
      intro _
      rcases s.eq_empty_or_nonempty with rfl | hs
      · exact ⟨a, by simp, by simp⟩
      · obtain ⟨j, hjs, hj⟩ := ih hs
        rcases ValuationRing.dvd_total (y a) (y j) with h | h
        · exact ⟨a, by simp, fun k hk => by
            rcases Finset.mem_insert.mp hk with rfl | hk'
            · exact dvd_rfl
            · exact h.trans (hj k hk')⟩
        · exact ⟨j, by simp [hjs], fun k hk => by
            rcases Finset.mem_insert.mp hk with rfl | hk'
            · exact h
            · exact hj k hk'⟩
  obtain ⟨j, _, hj⟩ := key Finset.univ ⟨Classical.arbitrary ι, Finset.mem_univ _⟩
  exact ⟨j, fun k => hj k (Finset.mem_univ _)⟩

/-- ★**1 つが他を全部割るなら、生成するイデアルはそれが生成する**。 -/
theorem span_range_eq_span_singleton_of_dvd {S : Type} [CommRing S] {ι : Type}
    (y : ι → S) (j : ι) (h : ∀ k, y j ∣ y k) :
    Ideal.span (Set.range y) = Ideal.span {y j} := by
  refine le_antisymm ?_ ?_
  · rw [Ideal.span_le]
    rintro _ ⟨k, rfl⟩
    exact Ideal.mem_span_singleton.mpr (h k)
  · rw [Ideal.span_le, Set.singleton_subset_iff]
    exact Ideal.subset_span ⟨j, rfl⟩

/-! ## ★★★★★★★★★★★★★★★★数体の整数環で -/

/-- ★★★★★★★★★★★★★★★★**有限素点 `Q` では座標の 1 つが `𝔞` を生成し、
`x_0/x_j` が `(𝓞_F)_Q` に入る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★これが `§9-940` の**有限素点の整合**の後ろ 2 つ
（`𝔞_Q = (x_j)` と `(s_0/s_j)(y)·x_j = x_0` の代数の部分）である。
★機構は `(𝓞_F)_Q` が離散付値環であること
（`IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain`）だけである。 -/
theorem exists_span_and_ratio_localization (F : Type) [Field F] [NumberField F]
    {N : ℕ} (x : Fin (N + 1) → 𝓞 F) (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal) :
    haveI := hQ.isPrime
    ∃ (j : Fin (N + 1)) (g : Localization.AtPrime Q),
      Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime Q)) (Ideal.span (Set.range x))
          = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)} ∧
        g * algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)
          = algebraMap (𝓞 F) (Localization.AtPrime Q) (x 0) := by
  haveI := hQ.isPrime
  haveI hQ0 : Q ≠ ⊥ := Ring.ne_bot_of_isMaximal_of_not_isField hQ (RingOfIntegers.not_isField F)
  haveI : IsDiscreteValuationRing (Localization.AtPrime Q) :=
    IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain (𝓞 F) hQ0 _
  obtain ⟨j, hj⟩ := exists_dvd_all_of_valuationRing
    (fun k => algebraMap (𝓞 F) (Localization.AtPrime Q) (x k))
  obtain ⟨g, hg⟩ := hj 0
  refine ⟨j, g, ?_, ?_⟩
  · rw [Ideal.map_span, ← Set.range_comp]
    exact span_range_eq_span_singleton_of_dvd _ j hj
  · rw [hg]
    ring

/-! ## ★出典の紐付け(`.src`) -/

def exists_dvd_all_of_valuationRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(付値環では有限族のどれか 1 つが他を全部割る)",
    sectionId := "genell-prop-1-4" }

def span_range_eq_span_singleton_of_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(1 つが他を全部割るなら生成するイデアルはそれが生成する)",
    sectionId := "genell-prop-1-4" }

def exists_span_and_ratio_localization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(有限素点では座標の 1 つが 𝔞 を生成する)",
    sectionId := "genell-prop-1-4" }

def exists_span_and_ratio_localization.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain((𝓞_F)_Q は DVR)"
      (.inMathlib "IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain") 2,
    .citation "[mathlib]" "ValuationRing.dvd_total(付値環では割り切りが全順序)"
      (.inMathlib "ValuationRing.dvd_total") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-940 の**有限素点の整合**のうち後ろ 2 つ" ++
       "(𝔞_Q = (x_j) と x_0/x_j ∈ (𝓞_F)_Q)は**代数だけで出る**。" ++
       "(𝓞_F)_Q は離散付値環で、付値環では割り切りが全順序だから、" ++
       "有限個の座標のうちどれか 1 つが他を全部割る" ++
       "——これが「v_Q(x_i) が最小の i を取る」の付値を経由しない形である") 5,
    .implicitStep
      ("★残るのは幾何の側の紐である——ここで得た j が §9-939 の局所チャートの添字と" ++
       "一致すること。機構は §9-913(ψ⁻¹(D₊(x_j)) = X_{s_j})と、" ++
       "座標が ψ から読まれていること(§9-941)の合わせ技である") 4 ]

end ABC3.Found.GenEll
