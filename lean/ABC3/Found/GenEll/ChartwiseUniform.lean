/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VerticalDescent

/-!
# 一様な `m` を**アフィン被覆から**取る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

## ★★★★★★★★★★★段 3b の後半 —— これで鎖が繋がる

`ResearchPaper/mathlib-gap.json` の `vertical-divisors-generate-pic-kernel` は
段 3b を 2 本に割っていた:

* **3b-1**（`VerticalDescent.lean`、済）: 層の水準の `n·D ≤ E` を点へ降ろす
* **3b-2**（★本ファイル）: その `n` を**アフィン被覆から**取る

★★本ファイルの到達点:

> **`X` が有限個のネーターなアフィンチャートで覆われ、
> 各チャートで `D` と `E` が `N` を反転して一致するなら、
> `ht_D ≈ ht_E`**（定数は `m·log N`）

## ★★★★★★機構は 4 段

| 段 | 内容 |
|---|---|
| 1 | 各チャートで `∃m_i, N^{m_i}·D_i ≤ E_i`（`IdealComparable.lean` の `exists_pow_mul_le_of_map_le`） |
| 2 | **有限**被覆なので `m := max m_i` が取れる |
| 3 | チャートごとの不等式を層の不等式へ（`sheaf_le_of_chartwise`） |
| 4 | 層の不等式を点へ（`VerticalDescent.lean`） |

★★★段 3 の要は `specIdealSheaf_comap_ideal_top` である——
`Spec ℤ` 上のイデアル層をアフィン `W` へ引き戻すと `J·Γ(W,⊤)` になる。
**`ℤ → Γ(W,⊤)` の環準同型は 1 つしかない**（`RingHom.ext_int`）から、
`W` にも `f` にも依らない。

## ★逸脱（明示）

★局所化のモデルを `Localization (Submonoid.powers N)` に固定した。
★★`IsLocalization` を満たす任意の環でも同じ結論になるので制限にならないが、
`∀ S [instances]` で量化すると `Fintype ι` 上での `choose` が煩雑になるため
具体のモデルを取った。

★★★ネーター性はチャートごとの仮定 `hnoeth` として受けている
（`X` が `ℤ` 上有限型ならすべて満たされる）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits NumberField

/-! ## ★★`≤` 版の貼り合わせ -/

/-- ★★**開埋め込み `U.ι` に沿った `comap` が `≤` なら、`U` での値も `≤`**。

★`ComapMul.lean` の `ideal_eq_of_comap_ι_eq`（等号版）の `≤` 版である。
★★機構も同じ——`appIso` の逆を `Iso.hom_inv_id` で圏の言葉のまま打ち消す。 -/
theorem ideal_le_of_comap_ι_le {X : Scheme.{0}} {A B : X.IdealSheafData} (U : X.affineOpens)
    (h : A.comap U.1.ι ≤ B.comap U.1.ι) : A.ideal U ≤ B.ideal U := by
  haveI : IsAffine U.1.toScheme := U.2
  have hUeq : (⟨U.1.ι ''ᵁ (⊤ : U.1.toScheme.Opens),
      (isAffineOpen_top U.1.toScheme).image_of_isOpenImmersion U.1.ι⟩ : X.affineOpens) = U := by
    apply Subtype.ext; simp
  have h2 : (A.comap U.1.ι).ideal ⟨⊤, isAffineOpen_top U.1.toScheme⟩
      ≤ (B.comap U.1.ι).ideal ⟨⊤, isAffineOpen_top U.1.toScheme⟩ :=
    Scheme.IdealSheafData.le_def.1 h _
  simp only [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion] at h2
  have hmono := Ideal.comap_mono (f := ((Scheme.Hom.appIso U.1.ι
    (⊤ : U.1.toScheme.Opens)).hom).hom) h2
  have hcomp : ∀ K : Ideal Γ(X, U.1.ι ''ᵁ (⊤ : U.1.toScheme.Opens)),
      Ideal.comap ((Scheme.Hom.appIso U.1.ι (⊤ : U.1.toScheme.Opens)).hom).hom
        (Ideal.comap ((Scheme.Hom.appIso U.1.ι (⊤ : U.1.toScheme.Opens)).inv).hom K) = K := by
    intro K
    rw [Ideal.comap_comap, ← CommRingCat.hom_comp, Iso.hom_inv_id]
    simp
  rw [hcomp, hcomp] at hmono
  rw [← hUeq]
  exact hmono

/-! ## ★★★`Spec ℤ` の層をアフィンチャートで読む -/

/-- ★★★★★**`Spec ℤ` 上のイデアル層をアフィン `W` へ引き戻すと `J·Γ(W,⊤)`**。

★★**`ℤ → Γ(W,⊤)` の環準同型は 1 つしかない**（`RingHom.ext_int`）ので、
結果は `W` にも `g` にも依らない。
★★★これが「`Spec ℤ` から来る分は一様」の、チャート版である。 -/
theorem specIdealSheaf_comap_ideal_top {W : Scheme.{0}} [IsAffine W]
    (g : W ⟶ Spec (CommRingCat.of ℤ)) (J : Ideal (CommRingCat.of ℤ)) :
    ((specIdealSheaf (CommRingCat.of ℤ) J).comap g).ideal ⟨⊤, isAffineOpen_top W⟩
      = J.map (Int.castRingHom Γ(W, ⊤)) := by
  rw [ideal_comap_eq_map_of_isAffine]
  show ((specIdealSheaf (CommRingCat.of ℤ) J).ideal
    ⟨⊤, isAffineOpen_top (Spec (CommRingCat.of ℤ))⟩).map (Scheme.Hom.appTop g).hom = _
  show (Scheme.IdealSheafData.equivOfIsAffine
      (Scheme.IdealSheafData.equivOfIsAffine.symm
        (J.map (Scheme.ΓSpecIso (CommRingCat.of ℤ)).inv.hom))).map (Scheme.Hom.appTop g).hom = _
  rw [OrderRingIso.apply_symm_apply, Ideal.map_map]
  congr 1
  exact RingHom.ext_int _ _

/-- ★**`specIdealSheaf` は単調**。 -/
theorem specIdealSheaf_mono {A : CommRingCat.{0}} {J J' : Ideal A} (h : J ≤ J') :
    specIdealSheaf A J ≤ specIdealSheaf A J' := by
  show Scheme.IdealSheafData.equivOfIsAffine.symm (J.map (Scheme.ΓSpecIso A).inv.hom)
    ≤ Scheme.IdealSheafData.equivOfIsAffine.symm (J'.map (Scheme.ΓSpecIso A).inv.hom)
  exact (map_le_map_iff Scheme.IdealSheafData.equivOfIsAffine.symm).2 (Ideal.map_mono h)

/-! ## ★★★★★★★★チャートごとの不等式を層の不等式へ -/

/-- ★★★★★★★★**チャートごとの `n`-不等式から層の `n`-不等式へ**。

★被覆の上で見れば足りる（mathlib の `le_of_iSup_eq_top`）。
★★各チャートでは `comap_mul`（`ComapMul.lean`）と
`specIdealSheaf_comap_ideal_top` で環の言葉に落ちる。 -/
theorem sheaf_le_of_chartwise {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    (D E : X.IdealSheafData) (n : ℤ)
    {ι : Type} (U : ι → X.affineOpens) (hcov : ⨆ i, ((U i : X.affineOpens) : X.Opens) = ⊤)
    (h : ∀ i, Ideal.span {((n : ℤ) : Γ((U i).1.toScheme, ⊤))}
        * (D.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩
      ≤ (E.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩) :
    (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap f * D ≤ E := by
  refine Scheme.IdealSheafData.le_of_iSup_eq_top U hcov (fun i => ?_)
  haveI : IsAffine (U i).1.toScheme := (U i).2
  apply ideal_le_of_comap_ι_le
  refine Scheme.IdealSheafData.le_of_isAffine ?_
  have hcm : ((specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap f * D).comap (U i).1.ι
      = (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap ((U i).1.ι ≫ f)
        * D.comap (U i).1.ι := by
    rw [ABC3.Found.GenEll.comap_mul, Scheme.IdealSheafData.comap_comp]
  rw [hcm]
  have hprod : ((specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap ((U i).1.ι ≫ f)
        * D.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top (U i).1.toScheme⟩
      = ((specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap
            ((U i).1.ι ≫ f)).ideal ⟨⊤, isAffineOpen_top (U i).1.toScheme⟩
        * (D.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top (U i).1.toScheme⟩ := by
    rw [Scheme.IdealSheafData.ideal_mul]; rfl
  rw [hprod, specIdealSheaf_comap_ideal_top, Ideal.map_span]
  have hcast : (Int.castRingHom Γ((U i).1.toScheme, ⊤)) n
      = ((n : ℤ) : Γ((U i).1.toScheme, ⊤)) := rfl
  simp only [Set.image_singleton, hcast]
  exact h i

/-! ## ★★★★★★★★★一様化 —— 被覆が**有限**だから最大が取れる -/

/-- ★★★★★★★★★**チャートごとの `m_i` を 1 つに揃える**。

★★**被覆が有限であることだけ**を使う（`Finset.univ.sup`）。
★これが「一様な `m`」の正体である。 -/
theorem exists_uniform_sheaf_le {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ))
    (D E : X.IdealSheafData) (N : ℤ)
    {ι : Type} [Fintype ι] [DecidableEq ι] (U : ι → X.affineOpens)
    (hcov : ⨆ i, ((U i : X.affineOpens) : X.Opens) = ⊤)
    (h : ∀ i, ∃ m : ℕ, Ideal.span {((N : Γ((U i).1.toScheme, ⊤))) ^ m}
        * (D.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩
      ≤ (E.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩) :
    ∃ m : ℕ, (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {N ^ m})).comap f * D ≤ E := by
  classical
  choose m hm using h
  refine ⟨Finset.univ.sup m, sheaf_le_of_chartwise f D E _ U hcov (fun i => ?_)⟩
  refine le_trans (Ideal.mul_mono_left ?_) (hm i)
  rw [Ideal.span_singleton_le_span_singleton]
  have hle : m i ≤ Finset.univ.sup m := Finset.le_sup (Finset.mem_univ i)
  refine ⟨((N : Γ((U i).1.toScheme, ⊤))) ^ (Finset.univ.sup m - m i), ?_⟩
  push_cast
  rw [← pow_add]
  congr 1
  omega

/-- ★★★★★★★★★★**「`N` を反転して一致する」から一様な `m` が出る**。

★各チャートは `IdealComparable.lean` の `exists_pow_mul_le_of_map_le`
（ネーター環で有限生成性だけを使う）が担当する。 -/
theorem exists_uniform_sheaf_le_of_localization {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (D E : X.IdealSheafData) (N : ℤ)
    {ι : Type} [Fintype ι] [DecidableEq ι] (U : ι → X.affineOpens)
    (hcov : ⨆ i, ((U i : X.affineOpens) : X.Opens) = ⊤)
    (hnoeth : ∀ i, IsNoetherianRing Γ((U i).1.toScheme, ⊤))
    (hagree : ∀ i,
      ((D.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers ((N : Γ((U i).1.toScheme, ⊤))))))
        ≤ ((E.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers ((N : Γ((U i).1.toScheme, ⊤))))))) :
    ∃ m : ℕ, (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {N ^ m})).comap f * D ≤ E := by
  refine exists_uniform_sheaf_le f D E N U hcov (fun i => ?_)
  haveI := hnoeth i
  exact exists_pow_mul_le_of_map_le ((N : Γ((U i).1.toScheme, ⊤)))
    (Localization (Submonoid.powers ((N : Γ((U i).1.toScheme, ⊤))))) _ _ (hagree i)

/-! ## ★★★★★★★★★★★到達点 —— 鎖が繋がった -/

/-- ★★★★★★★★★★★**チャートごとに `N` を反転して一致すれば、高さは BD-同値**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

> `X` が有限個のネーターなアフィンチャートで覆われ、
> 各チャートで `D` と `E` が `N` を反転して一致するなら **`ht_D ≈ ht_E`**。

★★★★★これが原文の `≈` の**完全な形式化**である。定数は `m·log N` で、
**点にも定義体にも依らない**。

★★機構の全体（本ファイル＋前 3 ファイル）:

1. 各チャートで `∃m_i, N^{m_i}·D_i ≤ E_i`（有限生成性のみ）
2. 有限被覆なので `m := max m_i`
3. チャートごとの不等式 → 層の不等式（`sheaf_le_of_chartwise`）
4. 層の不等式 → 点の不等式（`VerticalDescent.lean`。`Spec ℤ` が終対象）
5. 点の不等式 → 次数の差 ≤ `log n`（`VerticalBound.lean`。`absNorm` の単調性）
6. 次数の差 → 高さの BD-同値（`IdealComparable.lean`）

★★★★**交点数も Weil 因子類群も使っていない。** -/
theorem htArith_bdeq_of_chartwise_localization (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (D E : ArithCartier X)
    (N : ℕ) (hN : N ≠ 0)
    {ι : Type} [Fintype ι] [DecidableEq ι] (U : ι → X.affineOpens)
    (hcov : ⨆ i, ((U i : X.affineOpens) : X.Opens) = ⊤)
    (hnoeth : ∀ i, IsNoetherianRing Γ((U i).1.toScheme, ⊤))
    (hDE : ∀ i,
      ((D.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤))))))
        ≤ ((E.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤)))))))
    (hED : ∀ i,
      ((E.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤))))))
        ≤ ((D.divisor.comap (U i).1.ι).ideal ⟨⊤, isAffineOpen_top _⟩).map
          (algebraMap Γ((U i).1.toScheme, ⊤)
            (Localization (Submonoid.powers (((N : ℤ) : Γ((U i).1.toScheme, ⊤)))))))
    (hD0 : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0)
    (hE0 : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F E.divisor xF ≠ 0)
    (harc : ∀ xF : specRingOfIntegers F ⟶ X,
      (archADiv F D.green xF).sum (fun _ r => r)
        = (archADiv F E.green xF).sum (fun _ r => r)) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E xF) := by
  obtain ⟨m₁, h₁⟩ :=
    exists_uniform_sheaf_le_of_localization f D.divisor E.divisor (N : ℤ) U hcov hnoeth hDE
  obtain ⟨m₂, h₂⟩ :=
    exists_uniform_sheaf_le_of_localization f E.divisor D.divisor (N : ℤ) U hcov hnoeth hED
  have hsub : ∀ k : ℕ, k ≤ max m₁ m₂ →
      (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(N : ℤ) ^ max m₁ m₂})).comap f
        ≤ (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(N : ℤ) ^ k})).comap f := by
    intro k hk
    refine Scheme.IdealSheafData.comap_mono f (specIdealSheaf_mono ?_)
    rw [Ideal.span_singleton_le_span_singleton]
    exact ⟨(N : ℤ) ^ (max m₁ m₂ - k), by rw [← pow_add]; congr 1; omega⟩
  have hcast : ((N ^ max m₁ m₂ : ℕ) : ℤ) = (N : ℤ) ^ max m₁ m₂ := by push_cast; ring
  refine htArith_bdeq_of_sheaf_comparable F f D E (N ^ max m₁ m₂)
    (pow_ne_zero _ hN) ?_ ?_ hD0 hE0 harc
  · rw [hcast]; exact le_trans (mul_le_mul' (hsub m₁ (le_max_left _ _)) le_rfl) h₁
  · rw [hcast]; exact le_trans (mul_le_mul' (hsub m₂ (le_max_right _ _)) le_rfl) h₂

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (ii)` の証明中の 1 文
（`ht_{L⊗M} ≈ ht_L`）の**完全な形**である。命題 (ii) 全体ではない。 -/

def sheaf_le_of_chartwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——チャートごとの不等式を層へ)",
    sectionId := "genell-prop-1-4" }

def exists_uniform_sheaf_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——有限被覆で m を 1 つに揃える)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_chartwise_localization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——ht_{L⊗M} ≈ ht_L の完全な形)",
    sectionId := "genell-prop-1-4" }

def htArith_bdeq_of_chartwise_localization.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pow_mul_le_of_map_le(各チャートの m_i——有限生成性のみ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_mul_le_of_map_le") 6,
    .citation "[ABC3]" "span_natCast_mul_pullbackIdeal_le(層の仮定を点へ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.span_natCast_mul_pullbackIdeal_le") 6,
    .citation "[ABC3]" "comap_mul(2026-08-17 の平行セッション)"
      (.inProject "ABC3" "ABC3.Found.GenEll.comap_mul") 6,
    .citation "[mathlib]" "Scheme.IdealSheafData.le_of_iSup_eq_top(被覆の上で見れば足りる)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.le_of_iSup_eq_top") 6,
    .citation "[mathlib]" "RingHom.ext_int(ℤ からの環準同型は 1 つ)"
      (.inMathlib "RingHom.ext_int") 6,
    .implicitStep
      ("★★★★★原文は ht_{L⊗M} ≈ ht_L を 1 行で済ませている。" ++
       "★形式化では 6 段(チャートごとの m_i・有限被覆での最大・層への持ち上げ・" ++
       "点への引き戻し・次数の評価・高さへ)に分かれた。" ++
       "★★**交点数も Weil 因子類群も使っていない**——" ++
       "有限生成性と準コンパクト性だけである") 6,
    .implicitStep
      ("★逸脱: 局所化のモデルを Localization (Submonoid.powers N) に固定した。" ++
       "★★ネーター性はチャートごとの仮定 hnoeth として受けている" ++
       "(X が ℤ 上有限型ならすべて満たされる)") 6 ]

end ABC3.Found.GenEll
