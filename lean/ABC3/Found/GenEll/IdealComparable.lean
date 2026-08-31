/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VerticalBound
import ABC3.Found.GenEll.ConductorHeight

/-!
# [GenEll] `N` で比較できる 2 つのイデアル（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

## ★★★★★★★★★★仮定はここまで弱くできる

`VerticalTwist.lean` → `VerticalBound.lean` と弱めてきた仮定を、
**最も弱い形**にする:

> `N·I ⊆ J` かつ `N·J ⊆ I`  ⟹  `|deg I − deg J| ≤ log N`

★★これは「`I` と `J` は `N` を反転すれば同じ」という主張の**定量版**である。
★★★定数 `log N` は `I`・`J`・`F` のいずれにも依らない。

## ★★★★★★「`N` を反転すれば同じ」から `N^m` が出る

`exists_pow_mul_le_of_map_le`:

> `R` がネーター環で `I·R[1/N] ⊆ J·R[1/N]` なら、ある `m` で `N^m·I ⊆ J`。

★機構は初等的である——`I` は有限生成、各生成元 `a` について
`a/1 ∈ J·R[1/N]` から `N^{k_a}·a ∈ J` が出るので、`m = max k_a` を取る。

★★**これが「垂直な差は有界」の正体である**。交点数も類群も要らない。

## ★★★★★残っている段（明示）

★上の `m` は **`R` ごと**に決まる。点ごとに `𝓞_F` へ下ろすと `m` が点に依ってしまう。
★★**一様な `m` はスキームの側から来なければならない**——
`X` は準コンパクトなので**アフィン被覆は有限**であり、
各チャートの `m_i` の最大を取れば `X` 全体で 1 つの `m` が取れる。
★★★その段（イデアル層の水準での `N^m·I_D ⊆ I_E`）は本ファイルには入っていない。
`htArith_bdeq_of_ideal_comparable` はそれを**仮説として受けている**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField IsLocalization

/-! ## ★局所化で入るなら、`N` の冪を掛ければ入る -/

/-- ★★**`a/1 ∈ J·R[1/N]` なら `N^k·a ∈ J`**。

★機構は `IsLocalization.mem_map_algebraMap_iff` と `eq_iff_exists` の 2 本だけ。 -/
theorem exists_pow_mul_mem {R : Type} [CommRing R] (N : R) (S : Type) [CommRing S] [Algebra R S]
    [IsLocalization (Submonoid.powers N) S] (J : Ideal R) (a : R)
    (h : algebraMap R S a ∈ J.map (algebraMap R S)) :
    ∃ k : ℕ, N ^ k * a ∈ J := by
  rw [IsLocalization.mem_map_algebraMap_iff (Submonoid.powers N) S] at h
  obtain ⟨⟨⟨y, hy⟩, ⟨m, hm⟩⟩, hEq⟩ := h
  simp only at hEq
  rw [← map_mul] at hEq
  obtain ⟨⟨c, hc⟩, hcEq⟩ := (IsLocalization.eq_iff_exists (Submonoid.powers N) S).1 hEq
  obtain ⟨i, hi⟩ := hm
  obtain ⟨j, hj⟩ := hc
  simp only at hi hj hcEq
  refine ⟨j + i, ?_⟩
  have hkey : N ^ (j + i) * a = c * (a * m) := by
    rw [pow_add, ← hi, ← hj]; ring
  rw [hkey, hcEq]
  exact J.mul_mem_left _ hy

/-- ★★★★★★★★**ネーター環なら `m` を 1 つに揃えられる**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

`I·R[1/N] ⊆ J·R[1/N]` なら、ある `m` で `N^m·I ⊆ J`。

★★`I` が有限生成であることだけを使う。生成元ごとの `k_a` の**最大**を取る。
★★★これが「垂直な差は有界」の正体である——**交点数も類群も要らない**。 -/
theorem exists_pow_mul_le_of_map_le {R : Type} [CommRing R] [IsNoetherianRing R]
    (N : R) (S : Type) [CommRing S] [Algebra R S]
    [IsLocalization (Submonoid.powers N) S] (I J : Ideal R)
    (h : I.map (algebraMap R S) ≤ J.map (algebraMap R S)) :
    ∃ m : ℕ, Ideal.span {N ^ m} * I ≤ J := by
  classical
  obtain ⟨s, hs⟩ := (isNoetherian_def.mp inferInstance) I
  have hgen : ∀ a ∈ s, ∃ k : ℕ, N ^ k * a ∈ J := by
    intro a ha
    refine exists_pow_mul_mem N S J a (h ?_)
    exact Ideal.mem_map_of_mem _ (hs ▸ Ideal.subset_span ha)
  choose! k hk using hgen
  refine ⟨s.sup k, ?_⟩
  have hmem : ∀ x ∈ I, N ^ (s.sup k) * x ∈ J := by
    intro x hx
    rw [← hs] at hx
    refine Submodule.span_induction ?_ ?_ ?_ ?_ hx
    · intro a ha
      have hle : k a ≤ s.sup k := Finset.le_sup ha
      have hsplit : N ^ (s.sup k) * a = N ^ (s.sup k - k a) * (N ^ (k a) * a) := by
        rw [← mul_assoc, ← pow_add, Nat.sub_add_cancel hle]
      rw [hsplit]
      exact J.mul_mem_left _ (hk a ha)
    · simp
    · intro a b _ _ hA hB
      rw [mul_add]; exact J.add_mem hA hB
    · intro c a _ hA
      have hc : N ^ (s.sup k) * (c • a) = c * (N ^ (s.sup k) * a) := by
        rw [smul_eq_mul]; ring
      rw [hc]
      exact J.mul_mem_left _ hA
  refine Ideal.mul_le.2 fun r hr x hx => ?_
  obtain ⟨c, rfl⟩ := Ideal.mem_span_singleton'.1 hr
  rw [mul_assoc]
  exact J.mul_mem_left _ (hmem x hx)

section Degrees

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★次数の差は `log n` で抑えられる -/

/-- ★★**`n·I ⊆ J` なら `deg J ≤ deg I + log n`**。

★反単調性（`VerticalBound.lean`）と積の加法性（`DegMul.lean`）を繋ぐだけ。 -/
theorem degNormalized_idealADiv_le_add_log
    (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (n : ℕ) (hn : n ≠ 0)
    (h : Ideal.span {((n : ℕ) : 𝓞 F)} * I ≤ J) :
    degNormalized (idealADiv F J) ≤ degNormalized (idealADiv F I) + Real.log n := by
  have hne0 : ((n : ℕ) : 𝓞 F) ≠ 0 := Nat.cast_ne_zero.mpr hn
  have hspan0 : Ideal.span {((n : ℕ) : 𝓞 F)} ≠ 0 := by
    simpa [Ideal.span_singleton_eq_bot] using hne0
  have hmul0 : Ideal.span {((n : ℕ) : 𝓞 F)} * I ≠ 0 := mul_ne_zero hspan0 hI
  calc degNormalized (idealADiv F J)
      ≤ degNormalized (idealADiv F (Ideal.span {((n : ℕ) : 𝓞 F)} * I)) :=
        degNormalized_idealADiv_antitone F _ _ hmul0 h
    _ = degNormalized (idealADiv F (Ideal.span {((n : ℕ) : 𝓞 F)}))
          + degNormalized (idealADiv F I) :=
        degNormalized_idealADiv_mul F _ _ hspan0 hI
    _ = Real.log n + degNormalized (idealADiv F I) := by
        rw [degNormalized_idealADiv_span_natCast F n hn]
    _ = degNormalized (idealADiv F I) + Real.log n := by ring

/-- ★★★★★★★★**`n` で互いに比較できる 2 つのイデアルは次数が `log n` 以内**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★★これが「`N` を反転すれば同じ」の**定量版**である。
★★★定数は `I`・`J`・`F` のいずれにも依らない。 -/
theorem abs_degNormalized_idealADiv_sub_le
    (I J : Ideal (𝓞 F)) (hI : I ≠ 0) (hJ : J ≠ 0) (n : ℕ) (hn : n ≠ 0)
    (h1 : Ideal.span {((n : ℕ) : 𝓞 F)} * I ≤ J)
    (h2 : Ideal.span {((n : ℕ) : 𝓞 F)} * J ≤ I) :
    |degNormalized (idealADiv F I) - degNormalized (idealADiv F J)| ≤ Real.log n := by
  have hA := degNormalized_idealADiv_le_add_log F I J hI n hn h1
  have hB := degNormalized_idealADiv_le_add_log F J I hJ n hn h2
  rw [abs_le]
  constructor <;> linarith

/-! ## ★★★★★★★★★高さへの適用 -/

/-- ★**アルキメデス側が一致すれば、高さの差は有限素点側の次数の差**。 -/
theorem htArith_sub_eq_degNormalized_sub {X X' : Scheme.{0}}
    (D : ArithCartier X) (E : ArithCartier X')
    (xF : specRingOfIntegers F ⟶ X) (yF : specRingOfIntegers F ⟶ X')
    (harc : (archADiv F D.green xF).sum (fun _ r => r)
          = (archADiv F E.green yF).sum (fun _ r => r)) :
    htArith F D xF - htArith F E yF
      = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
        - degNormalized (idealADiv F (pullbackIdeal F E.divisor yF)) := by
  rw [htArith, htArith, degNormalized, degNormalized, degNormalized, degNormalized,
    deg_pullbackADiv, deg_pullbackADiv, harc]
  ring

/-- ★★★★★★★★★★**最も弱い仮定での高さの BD-同値**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

引き戻したイデアルが `n` で互いに比較でき、アルキメデス側が一致すれば、
**高さは BD-同値で定数は `log n`**。

| 版 | 仮定 |
|---|---|
| `VerticalTwist` | 差が `Spec ℤ` 上の因子の**底変換そのもの** |
| `VerticalBound` | 差のイデアルが**有理整数 `n` を含む** |
| **本定理** | 2 つのイデアルが **`n` で互いに比較できる** |

★★`n` が点に依らないことが本質である。
★★★その一様性は**スキームの準コンパクト性**から来る
（`exists_pow_mul_le_of_map_le` を有限個のアフィンチャートに当てる）——
★その段は本ファイルには入っていない。 -/
theorem htArith_bdeq_of_ideal_comparable {X X' : Scheme.{0}}
    (D : ArithCartier X) (E : ArithCartier X')
    (ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'))
    (n : ℕ) (hn : n ≠ 0)
    (hD0 : ∀ xF, pullbackIdeal F D.divisor xF ≠ 0)
    (hE0 : ∀ xF, pullbackIdeal F E.divisor (ePt xF) ≠ 0)
    (h1 : ∀ xF, Ideal.span {((n : ℕ) : 𝓞 F)} * pullbackIdeal F D.divisor xF
        ≤ pullbackIdeal F E.divisor (ePt xF))
    (h2 : ∀ xF, Ideal.span {((n : ℕ) : 𝓞 F)} * pullbackIdeal F E.divisor (ePt xF)
        ≤ pullbackIdeal F D.divisor xF)
    (harc : ∀ xF, (archADiv F D.green xF).sum (fun _ r => r)
          = (archADiv F E.green (ePt xF)).sum (fun _ r => r)) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E (ePt xF)) := by
  refine ⟨Real.log n, fun xF => ?_⟩
  show |htArith F D xF - htArith F E (ePt xF)| ≤ Real.log n
  rw [htArith_sub_eq_degNormalized_sub F D E xF (ePt xF) (harc xF)]
  exact abs_degNormalized_idealADiv_sub_le F _ _ (hD0 xF) (hE0 xF) n hn (h1 xF) (h2 xF)

end Degrees

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (ii)` の証明中の 1 段の**最弱形**である。命題全体ではない。 -/

def exists_pow_mul_le_of_map_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——N を反転して一致するなら N^m で比較できる)",
    sectionId := "genell-prop-1-4" }

def abs_degNormalized_idealADiv_sub_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——n で比較できるイデアルの次数差は log n 以内)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_ideal_comparable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——ht_{L⊗M} ≈ ht_L の最弱形)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_ideal_comparable.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degNormalized_idealADiv_antitone(次数はイデアルに反単調)"
      (.inProject "ABC3" "ABC3.Found.GenEll.degNormalized_idealADiv_antitone") 6,
    .citation "[ABC3]" "degNormalized_idealADiv_span_natCast(有理整数イデアルの次数は log n)"
      (.inProject "ABC3" "ABC3.Found.GenEll.degNormalized_idealADiv_span_natCast") 4,
    .citation "[mathlib]" "IsLocalization.mem_map_algebraMap_iff(局所化での像の判定)"
      (.inMathlib "IsLocalization.mem_map_algebraMap_iff") 6,
    .implicitStep
      ("★★★★★一様な n は**スキームの準コンパクト性**から来る" ++
       "——アフィン被覆が有限なので各チャートの m_i の最大を取れる。" ++
       "★exists_pow_mul_le_of_map_le が各チャートを担当する。" ++
       "★★その組み立て(イデアル層の水準)は本ファイルには入っていない") 6,
    .implicitStep
      ("★★アルキメデス側が一致することは、ℚ-同型から ℂ-点が対応することの帰結である。" ++
       "★本ファイルはそれを仮説 harc として受けている") 6 ]

end ABC3.Found.GenEll
