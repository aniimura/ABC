import ABC3.Found.GaloisRep.SupportSum
import ABC3.Found.GaloisRep.ClassSum

/-!
# Galois (G5) 第 174 ブロック —— **★★★★★★★★★D2 が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★`J` は単項である

第 165 の `rootUnit`(明示的な `n` 乗根イデアル)が**単項**であることを示す。
★これで (G5) 層 3 の `sorry` が消える。

    [J] = toClass( Σ_{v} e_v · Q_v )        (第 166)
        = toClass( α·(Σ_{[n]⁻¹(P)} Q) + β·(Σ_{T ∈ E[n]∖{O}} T) )
        = toClass( α·0 + β·0 ) = 0

### ★★★★★★台の二分法

`count_v ≠ 0` なら **`n·Q_v = P` か `n·Q_v = 0`** のどちらか(`nsmul_pointOf_dichotomy`)。
★前者は第 164(`count ≠ 0 ⟺ I_P = P'`)、後者は場合 B である。
★★第 171 で `hnn` が `n·Q_v ≠ 0` から自動になったので、この二分法が点の言葉で書ける。

### ★★★★★係数の一致は要らない

ファイバー上で `e_v = α`、`E[n]∖{O}` 上で `e_v = β`(第 171 の `count_eq_of_sub_torsion`)。
★**α = β を示す必要はない**——2 つの総和がそれぞれ独立に 0 だからである。
★★つまり **`[n]` の不分岐性(`e = 1`)を一度も使っていない**。

### ★★★★素点の総和を点の総和に移す

第 172 の往復(`pointOf`・`placeOf`)で `Finset` を移す:

    s := t.attach.image (placeOf ·)   ⟹   s.image pointOf = t

★`Finset.sum_image`(`pointOf` は単射)で `Σ_{v ∈ s}` が `Σ_{Q ∈ t}` になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `epsPoint` / `epsPoint_pointOf` | ★★点に付いた指数 |
| `pointOf_injective` / `placeOf_pointOf'` | ★★往復の帰結 |
| `image_pointOf_placeOf` | ★★★点の `Finset` ↔ 素点の `Finset` |
| `sum_places_eq_sum_points` | ★★★★★★**素点の総和を点の総和に** |
| `eps_eq_of_sub_torsion` | ★★★★★★★点の言葉での `e_v` の一定性 |
| `nsmul_pointOf_dichotomy` | ★★★★★★**台の二分法** |
| `sum_eps_split_zero` | ★★★★★★★★台の総和は 0 |
| `exists_spanSingleton_rootUnit_mulByN` | ★★★★★★★★★**D2 —— `J` は単項** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep
open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [IsAlgClosed F] [DecidableEq F]
  (W : WeierstrassCurve.Affine F) [W.IsElliptic]
  [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★★点に付いた指数 -/

/-- ★★点に付いた指数(0 でない点にはその素点の指数、無限遠点には 0)。 -/
noncomputable def epsPoint (I : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) (n : ℕ)
    (Q : W.Point) : ℤ :=
  if h : Q = 0 then 0
  else FractionalIdeal.count W.FunctionField (placeOf W Q h) I / (n : ℤ)

theorem epsPoint_pointOf (h2 : IsUnit (2 : F))
    (I : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) (n : ℕ)
    (v : HeightOneSpectrum W.CoordinateRing) :
    epsPoint W I n (pointOf W h2 v) = FractionalIdeal.count W.FunctionField v I / (n : ℤ) := by
  rw [epsPoint, dif_neg (pointOf_ne_zero W h2 v), placeOf_pointOf W h2 v]

theorem pointOf_injective (h2 : IsUnit (2 : F)) : Function.Injective (pointOf W h2) := by
  intro a b hab
  rw [← placeOf_pointOf W h2 a, ← placeOf_pointOf W h2 b]
  simp only [hab]

/-- ★★証明無関係版(`placeOf` の値は `≠ 0` の証明に依らない)。 -/
theorem placeOf_pointOf' (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (hne : pointOf W h2 v ≠ 0) : placeOf W (pointOf W h2 v) hne = v :=
  HeightOneSpectrum.ext (pointAsIdeal W h2 v).symm

/-! ## ★★★点の `Finset` と素点の `Finset` -/

variable [DecidableEq W.Point] [DecidableEq (HeightOneSpectrum W.CoordinateRing)]

/-- ★★★点の `Finset` から素点の `Finset` を作り、戻すと元に戻る。 -/
theorem image_pointOf_placeOf (h2 : IsUnit (2 : F)) (t : Finset W.Point)
    (hne : ∀ Q ∈ t, Q ≠ 0) :
    (t.attach.image (fun q : {Q // Q ∈ t} => placeOf W q.1 (hne q.1 q.2))).image
      (pointOf W h2) = t := by
  rw [Finset.image_image]
  have hcongr : ∀ q ∈ t.attach,
      ((pointOf W h2) ∘ (fun q : {Q // Q ∈ t} => placeOf W q.1 (hne q.1 q.2))) q = q.1 :=
    fun q _ => pointOf_placeOf W h2 q.1 (hne q.1 q.2)
  rw [Finset.image_congr hcongr, Finset.attach_image_val]

/-- ★★★★★★**素点の総和を点の総和に移す**。 -/
theorem sum_places_eq_sum_points (h2 : IsUnit (2 : F))
    (I : FractionalIdeal W.CoordinateRing⁰ W.FunctionField) (n : ℕ)
    (t : Finset W.Point) (hne : ∀ Q ∈ t, Q ≠ 0) :
    ∑ v ∈ t.attach.image (fun q : {Q // Q ∈ t} => placeOf W q.1 (hne q.1 q.2)),
        (FractionalIdeal.count W.FunctionField v I / (n : ℤ)) • pointOf W h2 v
      = ∑ Q ∈ t, epsPoint W I n Q • Q := by
  set s := t.attach.image (fun q : {Q // Q ∈ t} => placeOf W q.1 (hne q.1 q.2)) with hs
  have hstep : ∀ v ∈ s,
      (FractionalIdeal.count W.FunctionField v I / (n : ℤ)) • pointOf W h2 v
        = epsPoint W I n (pointOf W h2 v) • pointOf W h2 v :=
    fun v _ => by rw [epsPoint_pointOf W h2 I n v]
  have himg : ∑ Q ∈ s.image (pointOf W h2), epsPoint W I n Q • Q
      = ∑ x ∈ s, epsPoint W I n (pointOf W h2 x) • pointOf W h2 x :=
    Finset.sum_image (fun a _ b _ hab => pointOf_injective W h2 hab)
  rw [Finset.sum_congr rfl hstep, ← himg, hs, image_pointOf_placeOf W h2 t hne]

/-! ## ★★★★★★台の二分法と一定性 -/

/-- ★★★★★★**台の二分法**——`count ≠ 0` なら `n·Q_v = P` か `n·Q_v = 0`。 -/
theorem nsmul_pointOf_dichotomy (h2 : IsUnit (2 : F)) {n : ℕ} (hn : 1 ≤ n)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    {x y : F} (hP : W.Nonsingular x y) (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (v : HeightOneSpectrum W.CoordinateRing)
    (hcount : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) ≠ 0) :
    n • pointOf W h2 v = Point.some x y hP ∨ n • pointOf W h2 v = 0 := by
  by_cases hz : n • pointOf W h2 v = 0
  · exact Or.inr hz
  · refine Or.inl ?_
    obtain ⟨xR, yR, hR, hRm⟩ := exists_some_of_ne_zero W _ hz
    have hnn := valuation_mulByN_le_one W v (pointEq W h2 v) (pointAsIdeal W h2 v) h2 n μ hμF
      hns hμx hμy hμP hz
    have hkey := (count_ne_zero_iff_nsmul W v (pointEq W h2 v) (pointAsIdeal W h2 v) h2 hn μ
      hμinj hμF hns hμx hμy hμP hnn hP fP hfP hR hRm).1 hcount
    exact hRm.trans hkey.symm

variable [Infinite F]

/-- ★★★★★★★点の言葉での `e_v` の一定性。 -/
theorem eps_eq_of_sub_torsion (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (r : W.CoordinateRing) (hz : μ r ≠ 0)
    (Q Q' : W.Point) (hQ : Q ≠ 0) (hQ' : Q' ≠ 0) (hdiff : n • (Q' - Q) = 0) :
    epsPoint W (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r)) n Q
      = epsPoint W (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r)) n Q' := by
  obtain ⟨a, b, hb, rfl⟩ := exists_some_of_ne_zero W Q hQ
  obtain ⟨a', b', hb', rfl⟩ := exists_some_of_ne_zero W Q' hQ'
  rw [epsPoint, epsPoint, dif_neg hQ, dif_neg hQ']
  congr 1
  exact count_eq_of_sub_torsion W h2 h4 _ _ hb.1 (placeOf_some W hb hQ) hb'.1
    (placeOf_some W hb' hQ') n μ hμF hns hμx hμy hμP hdiff r hz

/-- ★★★★★★★★**台の総和は 0**——係数の一致は要らない。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★ファイバーと `E[n]∖{O}` の総和がそれぞれ独立に 0 だからである。 -/
theorem sum_eps_split_zero (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)]
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    {P : W.Point} (hPne : P ≠ 0)
    (s₁ : Finset W.Point)
    (hs₁sub : ∀ Q ∈ s₁, n • Q = P) (hs₁sum : ∑ Q ∈ s₁, Q = 0)
    (r : W.CoordinateRing) (hz : μ r ≠ 0) :
    ∑ Q ∈ s₁ ∪ ((Finset.univ : Finset (torsionPoints W n)).image
        (fun T : torsionPoints W n => (T : W.Point))).erase 0,
      epsPoint W (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r)) n Q • Q = 0 := by
  set s₂ := ((Finset.univ : Finset (torsionPoints W n)).image
    (fun T : torsionPoints W n => (T : W.Point))).erase 0 with hs₂def
  set ε := epsPoint W (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ r)) n with hεdef
  have hs₂sub : ∀ Q ∈ s₂, n • Q = 0 ∧ Q ≠ 0 := by
    intro Q hQ
    refine ⟨?_, Finset.ne_of_mem_erase hQ⟩
    obtain ⟨T, -, rfl⟩ := Finset.mem_image.1 (Finset.mem_of_mem_erase hQ)
    exact T.2
  have hconst : ∀ (u : Finset W.Point), (∀ Q ∈ u, Q ≠ 0) →
      (∀ Q ∈ u, ∀ Q' ∈ u, n • (Q' - Q) = 0) → ∃ α : ℤ, ∀ Q ∈ u, ε Q = α := by
    intro u hune hutor
    rcases u.eq_empty_or_nonempty with rfl | ⟨Q₀, hQ₀⟩
    · exact ⟨0, by simp⟩
    · refine ⟨ε Q₀, fun Q hQ => ?_⟩
      exact (eps_eq_of_sub_torsion W h2 h4 n μ hμF hns hμx hμy hμP r hz Q₀ Q
        (hune Q₀ hQ₀) (hune Q hQ) (hutor Q₀ hQ₀ Q hQ)).symm
  obtain ⟨α, hα⟩ := hconst s₁ (fun Q hQ => by
      intro h0; exact hPne (by rw [← hs₁sub Q hQ, h0, smul_zero]))
    (fun Q hQ Q' hQ' => by rw [smul_sub, hs₁sub Q' hQ', hs₁sub Q hQ, sub_self])
  obtain ⟨β, hβ⟩ := hconst s₂ (fun Q hQ => (hs₂sub Q hQ).2)
    (fun Q hQ Q' hQ' => by rw [smul_sub, (hs₂sub Q' hQ').1, (hs₂sub Q hQ).1, sub_self])
  refine sum_split_eq_zero (s₁ ∪ s₂) s₁ s₂ ε α β (Finset.Subset.refl _) ?_
    (fun Q hQ hQs => absurd hQ hQs) hα hβ hs₁sum (sum_torsion_erase_zero W n hn hchar)
  refine Finset.disjoint_left.2 (fun Q hQ1 hQ2 => ?_)
  exact hPne (by rw [← hs₁sub Q hQ1, (hs₂sub Q hQ2).1])

/-! ## ★★★★★★★★★D2 -/

/-- ★★★★★★★★★**D2 —— `J` は単項である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★台を `[n]⁻¹(P)` と `E[n]∖{O}` に分け、それぞれの点の総和が 0 であることから
`[J] = toClass(0) = 0`。★★分岐指数の値は一度も要らない。 -/
theorem exists_spanSingleton_rootUnit_mulByN (h2 : IsUnit (2 : F)) (h4 : (4 : F) ≠ 0)
    (n : ℕ) (hn : 1 ≤ n) (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    [Fintype (torsionPoints W n)]
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    {x y : F} (hP : W.Nonsingular x y) (hPt : n • Point.some x y hP = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) :
    ∃ g : W.FunctionField,
      ((rootUnit (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) n :
          (FractionalIdeal W.CoordinateRing⁰ W.FunctionField)ˣ)
        : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
      = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g := by
  set I := FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) with hIdef
  set PP := Point.some x y hP with hPPdef
  have hPne : PP ≠ 0 := by rw [hPPdef]; simp
  have hμfP : μ fP ≠ 0 := fun h0 =>
    (fP_ne_zero W n fP hfP) (hμinj (by rw [h0, map_zero]))
  set s₂ := ((Finset.univ : Finset (torsionPoints W n)).image
    (fun T : torsionPoints W n => (T : W.Point))).erase 0 with hs₂def
  obtain ⟨s₁, hs₁sub, hs₁sum, hs₁mem⟩ :
      ∃ s₁ : Finset W.Point, (∀ Q ∈ s₁, n • Q = PP) ∧ (∑ Q ∈ s₁, Q = 0)
        ∧ (∀ Q : W.Point, n • Q = PP → Q ∈ s₁) := by
    by_cases hex : ∃ Q₀ : W.Point, n • Q₀ = PP
    · obtain ⟨Q₀, hQ₀⟩ := hex
      refine ⟨(Finset.univ : Finset (torsionPoints W n)).image
        (fun T : torsionPoints W n => Q₀ + (T : W.Point)), ?_, ?_, ?_⟩
      · rintro Q hQ
        obtain ⟨T, -, rfl⟩ := Finset.mem_image.1 hQ
        rw [smul_add, hQ₀, T.2, add_zero]
      · exact sum_fiber_zero W n hn hchar Q₀ hQ₀ (hPPdef ▸ hPt)
      · intro Q hQ
        refine Finset.mem_image.2 ⟨⟨Q - Q₀, ?_⟩, Finset.mem_univ _, by simp⟩
        show n • (Q - Q₀) = 0
        rw [smul_sub, hQ, hQ₀, sub_self]
    · exact ⟨∅, by simp, by simp, fun Q hQ => absurd ⟨Q, hQ⟩ hex⟩
  set t := s₁ ∪ s₂ with htdef
  have hne : ∀ Q ∈ t, Q ≠ 0 := by
    intro Q hQ
    rcases Finset.mem_union.1 hQ with hh | hh
    · intro h0; exact hPne (by rw [← hs₁sub Q hh, h0, smul_zero])
    · exact Finset.ne_of_mem_erase hh
  set s := t.attach.image (fun q : {Q // Q ∈ t} => placeOf W q.1 (hne q.1 q.2)) with hsdef
  have hs : ∀ v : HeightOneSpectrum W.CoordinateRing,
      FractionalIdeal.count W.FunctionField v I ≠ 0 → v ∈ s := by
    intro v hcv
    have hd := nsmul_pointOf_dichotomy W h2 hn μ hμinj hμF hns hμx hμy hμP hP fP hfP v hcv
    have hmem : pointOf W h2 v ∈ t := by
      rcases hd with hh | hh
      · exact Finset.mem_union_left _ (hs₁mem _ hh)
      · refine Finset.mem_union_right _ (Finset.mem_erase.2 ⟨pointOf_ne_zero W h2 v, ?_⟩)
        exact Finset.mem_image.2 ⟨⟨pointOf W h2 v, hh⟩, Finset.mem_univ _, rfl⟩
    have himg : placeOf W (pointOf W h2 v) (hne _ hmem) ∈ s :=
      Finset.mem_image.2 ⟨⟨pointOf W h2 v, hmem⟩, Finset.mem_attach _ _, rfl⟩
    rwa [placeOf_pointOf' W h2 v (hne _ hmem)] at himg
  have hsum : ∑ v ∈ s, (FractionalIdeal.count W.FunctionField v I / (n : ℤ)) • pointOf W h2 v
      = 0 := by
    rw [hsdef, sum_places_eq_sum_points W h2 I n t hne, htdef]
    exact sum_eps_split_zero W h2 h4 n hn hchar μ hμF hns hμx hμy hμP hPne s₁ hs₁sub hs₁sum
      fP hμfP
  exact exists_spanSingleton_rootUnit W (pointOf W h2) (classGroup_primeUnit_pointOf W h2) I n s
    hs hsum

/-! ## ★出典の紐付け(`.src`) -/

def nsmul_pointOf_dichotomy.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——因子の台の二分法)",
    sectionId := "genell-thm-3-8" }

def exists_spanSingleton_rootUnit_mulByN.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——D2、n 乗根イデアルが単項であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
