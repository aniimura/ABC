/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Remark151Chain
import ABC3.Found.GenEll.ComapMul
import ABC3.Found.GenEll.IdealComparable

/-!
# 層の水準の `n·D ≤ E` を**点へ降ろす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

## ★★★★★★★★★★段 3b の前半

`ResearchPaper/mathlib-gap.json` の `vertical-divisors-generate-pic-kernel` の
段 3b は 2 本に分かれていた:

1. **`Spec ℤ` 上の `N` のイデアル層を点へ引き戻すと `span{N}` になること** ← ★本ファイル
2. 各チャートで `Γ(X_N ∩ U) = Γ(X,U)[1/N]` と、有限被覆で `m` の最大を取ること

★本ファイルは 1 を閉じ、**層の水準の仮定を点の水準の仮定に翻訳する**:

    `(Spec ℤ の n のイデアル層)^* · D ≤ E`
      ⟹  `∀ x,  span{n}·(x^*D) ≤ x^*E`

★★これで `IdealComparable.lean` の `htArith_bdeq_of_ideal_comparable` が
**層の水準の仮定から使える**ようになる。

## ★★★★★★機構は 3 本

| 補題 | 中身 | 出どころ |
|---|---|---|
| `pullbackIdealOf_mul` | 引き戻しは積を保つ | `ComapMul.lean` の `comap_mul`（平行セッション） |
| `pullbackIdealOf_mono` | 引き戻しは単調 | mathlib の `comap_mono` |
| `pullbackIdealOf_comap_intSpan` | `Spec ℤ` から来た層は `span{n}` に落ちる | `Spec ℤ` が**終対象**であること |

★★★3 本目が本質である——`x ≫ f` は `Spec ℤ` への**唯一の**射なので、
`X` にも点にも依らずに `ℤ → B` の構造射になる。

## ★残っている段（明示）

★段 3b の 2 本目（一様な `m` をアフィン被覆から取る段）は本ファイルには無い。
★★`IdealComparable.lean` の `exists_pow_mul_le_of_map_le` が各チャートを担当し、
`X` が準コンパクトなので有限個の最大を取ればよい、というところまでは判っている。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits NumberField

/-! ## ★引き戻しは積を保ち、単調である -/

/-- ★★★**点に沿った引き戻しは積を保つ**。

★`ComapMul.lean` の `comap_mul`（イデアル層の水準、平行セッションが 2026-08-17 に取得）を
`equivOfIsAffine`（環同型）と `Ideal.comap_symm` で環の水準へ落とすだけ。 -/
theorem pullbackIdealOf_mul {B : CommRingCat.{0}} {X : Scheme.{0}}
    (D E : X.IdealSheafData) (x : Spec B ⟶ X) :
    pullbackIdealOf B (D * E) x = pullbackIdealOf B D x * pullbackIdealOf B E x := by
  have keyB : ∀ J : Ideal Γ(Spec B, ⊤),
      Ideal.comap (Scheme.ΓSpecIso B).inv.hom J = J.map (Scheme.ΓSpecIso B).hom.hom :=
    fun J => Ideal.comap_symm (Scheme.ΓSpecIso B).commRingCatIsoToRingEquiv
  show Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine ((D * E).comap x)) = _
  rw [ABC3.Found.GenEll.comap_mul, map_mul, keyB]
  show _ = Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (D.comap x))
    * Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (E.comap x))
  rw [keyB, keyB, Ideal.map_mul]

/-- ★★**点に沿った引き戻しは単調**。 -/
theorem pullbackIdealOf_mono {B : CommRingCat.{0}} {X : Scheme.{0}}
    {D E : X.IdealSheafData} (h : D ≤ E) (x : Spec B ⟶ X) :
    pullbackIdealOf B D x ≤ pullbackIdealOf B E x := by
  show Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine (D.comap x)) ≤ _
  exact Ideal.comap_mono
    ((map_le_map_iff Scheme.IdealSheafData.equivOfIsAffine).2
      (Scheme.IdealSheafData.comap_mono x h))

/-! ## ★★★★アフィンスキーム上のイデアル層 -/

/-- ★★**イデアル `J ⊆ A` が定める `Spec A` 上のイデアル層**。

★`equivOfIsAffine`（`X.IdealSheafData ≃+*o Ideal Γ(X,⊤)`）の逆像として作る。 -/
noncomputable def specIdealSheaf (A : CommRingCat.{0}) (J : Ideal A) : (Spec A).IdealSheafData :=
  Scheme.IdealSheafData.equivOfIsAffine.symm (J.map (Scheme.ΓSpecIso A).inv.hom)

/-- ★**恒等射に沿った引き戻しは元のイデアル**（構成の検算）。 -/
theorem pullbackIdealOf_specIdealSheaf_id (A : CommRingCat.{0}) (J : Ideal A) :
    pullbackIdealOf A (specIdealSheaf A J) (𝟙 (Spec A)) = J := by
  show Ideal.comap _ (Scheme.IdealSheafData.equivOfIsAffine
    ((specIdealSheaf A J).comap (𝟙 (Spec A)))) = _
  rw [Scheme.IdealSheafData.comap_id]
  show Ideal.comap (Scheme.ΓSpecIso A).inv.hom
    (Scheme.IdealSheafData.equivOfIsAffine
      (Scheme.IdealSheafData.equivOfIsAffine.symm (J.map (Scheme.ΓSpecIso A).inv.hom))) = _
  rw [OrderRingIso.apply_symm_apply]
  exact Ideal.comap_map_of_bijective _
    (ConcreteCategory.bijective_of_isIso (Scheme.ΓSpecIso A).inv)

/-- ★★★★★★**`Spec.map φ` に沿った引き戻しは `Ideal.map φ`**。 -/
theorem pullbackIdealOf_specIdealSheaf {A B : CommRingCat.{0}} (J : Ideal A) (φ : A ⟶ B) :
    pullbackIdealOf B (specIdealSheaf A J) (Spec.map φ) = J.map φ.hom := by
  have h := pullbackIdealOf_specMap (B := A) (B' := B) (specIdealSheaf A J) (𝟙 (Spec A)) φ
  rw [Category.comp_id] at h
  rw [h, pullbackIdealOf_specIdealSheaf_id]

/-! ## ★★★★★★★★`Spec ℤ` から来た層は点で `span{n}` になる -/

/-- ★★★★★★★★**`Spec ℤ` 上のイデアル層を任意の点へ引き戻すと、係数の拡大イデアル**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★★機構は **`Spec ℤ` が終対象**であることだけである——
`x ≫ f` は `Spec ℤ` への**唯一の**射なので、`X` にも点にも依らない。
★★★これが「`Spec ℤ` から来る分は点に依らない」の形式化である。 -/
theorem pullbackIdealOf_comap_specIdealSheaf {B : CommRingCat.{0}} {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (J : Ideal (CommRingCat.of ℤ)) (x : Spec B ⟶ X) :
    pullbackIdealOf B ((specIdealSheaf (CommRingCat.of ℤ) J).comap f) x
      = J.map (CommRingCat.ofHom (Int.castRingHom B)).hom := by
  rw [← pullbackIdealOf_target]
  have hterm : x ≫ f = Spec.map (CommRingCat.ofHom (Int.castRingHom B)) :=
    specZIsTerminal.hom_ext _ _
  rw [hterm, pullbackIdealOf_specIdealSheaf]

/-- ★★★★★★★★★**有理整数 `n` の場合** —— `span {n}` になる。 -/
theorem pullbackIdealOf_comap_intSpan {B : CommRingCat.{0}} {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℤ) (x : Spec B ⟶ X) :
    pullbackIdealOf B
        ((specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap f) x
      = Ideal.span {((n : ℤ) : B)} := by
  rw [pullbackIdealOf_comap_specIdealSheaf, Ideal.map_span]
  congr 1
  simp

/-! ## ★★★★★★★★★★到達点 —— 層の仮定が点の仮定になる -/

/-- ★★★★★★★★★★**層の水準の `n·D ≤ E` は、点の水準へそのまま降りる**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

★★仮定は**層 1 つ**（`X` 上の不等式）で、結論は**すべての点**についての不等式である。
★★★定数 `n` が点に依らないのは、`Spec ℤ` が終対象だからである。 -/
theorem span_mul_pullbackIdealOf_le {B : CommRingCat.{0}} {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℤ) (D E : X.IdealSheafData)
    (h : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {n})).comap f * D ≤ E)
    (x : Spec B ⟶ X) :
    Ideal.span {((n : ℤ) : B)} * pullbackIdealOf B D x ≤ pullbackIdealOf B E x := by
  have h1 := pullbackIdealOf_mono h x
  rwa [pullbackIdealOf_mul, pullbackIdealOf_comap_intSpan] at h1

/-- ★★★★★★★★★★**数体の場合** ——
`IdealComparable.lean` の `htArith_bdeq_of_ideal_comparable` が要求する形そのもの。

★★これで「高さの BD-同値」の仮定が**層 1 つ**になった。 -/
theorem span_natCast_mul_pullbackIdeal_le (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} (f : X ⟶ Spec (CommRingCat.of ℤ)) (n : ℕ) (D E : X.IdealSheafData)
    (h : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(n : ℤ)})).comap f * D ≤ E)
    (xF : specRingOfIntegers F ⟶ X) :
    Ideal.span {((n : ℕ) : 𝓞 F)} * pullbackIdeal F D xF ≤ pullbackIdeal F E xF := by
  have h1 := span_mul_pullbackIdealOf_le (B := CommRingCat.of (𝓞 F)) f (n : ℤ) D E h xF
  have hcast : (((n : ℤ) : CommRingCat.of (𝓞 F))) = ((n : ℕ) : 𝓞 F) := by push_cast; ring
  rw [hcast] at h1
  exact h1

/-! ## ★★★★★★★★★★★高さへ —— 仮定が**層 1 つ**になった -/

/-- ★★★★★★★★★★★**層の水準で `n`-比較できる 2 つの算術因子の高さは BD-同値**。

原文 (GenEll p.6):
> Now observe that if M is an arithmetic line bundle that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z), then

`IdealComparable.lean` の `htArith_bdeq_of_ideal_comparable` は
**すべての点について**イデアルが `n`-比較できることを仮定していた。
★★本定理はそれを**層の不等式 2 本**に置き換える:

    `(Spec ℤ の n)^*·D ≤ E`  かつ  `(Spec ℤ の n)^*·E ≤ D`
      ⟹  `ht_D ≈ ht_E`（定数 `log n`）

★★★これが**幾何が実際に与える形**である
——`D` と `E` が `ℤ[1/n]` の上で一致すれば、両向きの不等式が
（`X` が準コンパクト・ネーターなら）ある `n` で成り立つ。

★★★★仮定 `harc`（アルキメデス側が一致）は、
`D` と `E` が同じ `ℂ`-点の上で同じ計量を持つことに対応する。 -/
theorem htArith_bdeq_of_sheaf_comparable (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (D E : ArithCartier X) (n : ℕ) (hn : n ≠ 0)
    (h1 : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(n : ℤ)})).comap f * D.divisor
        ≤ E.divisor)
    (h2 : (specIdealSheaf (CommRingCat.of ℤ) (Ideal.span {(n : ℤ)})).comap f * E.divisor
        ≤ D.divisor)
    (hD0 : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0)
    (hE0 : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F E.divisor xF ≠ 0)
    (harc : ∀ xF : specRingOfIntegers F ⟶ X,
      (archADiv F D.green xF).sum (fun _ r => r)
        = (archADiv F E.green xF).sum (fun _ r => r)) :
    BDeq (fun xF => htArith F D xF) (fun xF => htArith F E xF) :=
  htArith_bdeq_of_ideal_comparable F D E (fun xF => xF) n hn hD0 hE0
    (fun xF => span_natCast_mul_pullbackIdeal_le F f n D.divisor E.divisor h1 xF)
    (fun xF => span_natCast_mul_pullbackIdeal_le F f n E.divisor D.divisor h2 xF)
    harc

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (ii)` の証明中の段（垂直な差）を、
**層の水準から点の水準へ降ろす配管**である。 -/

def htArith_bdeq_of_sheaf_comparable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——層の水準で n-比較できれば高さは BD-同値)",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(引き戻しが積を保つこと——点の水準)",
    sectionId := "genell-prop-1-4" }

def pullbackIdealOf_comap_intSpan.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——Spec ℤ から来た層は点で span{n} になる)",
    sectionId := "genell-prop-1-4" }

def span_natCast_mul_pullbackIdeal_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(証明中の段——層の仮定を点の仮定へ降ろす)",
    sectionId := "genell-prop-1-4" }

def span_natCast_mul_pullbackIdeal_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "comap_mul(イデアル層の水準。2026-08-17 の平行セッション)"
      (.inProject "ABC3" "ABC3.Found.GenEll.comap_mul") 6,
    .citation "[mathlib]" "Scheme.IdealSheafData.comap_mono"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap_mono") 6,
    .citation "[mathlib]" "specZIsTerminal(Spec ℤ は終対象)"
      (.inMathlib "AlgebraicGeometry.specZIsTerminal") 6,
    .implicitStep
      ("★★★定数 n が点にも定義体にも依らないのは、" ++
       "x ≫ f が Spec ℤ への**唯一の**射だからである。" ++
       "★これが原文の「Spec(Z) から来る分」の形式化である") 6,
    .implicitStep
      ("★★残る段(段 3b の 2 本目): 一様な m をアフィン被覆から取ること。" ++
       "★各チャートは IdealComparable.lean の exists_pow_mul_le_of_map_le が担当し、" ++
       "X が準コンパクトなので有限個の最大を取ればよい") 6 ]

end ABC3.Found.GenEll
