/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.Cartier
import ABC3.Found.Divisor.CartierFrobenioid
import ABC3.Found.Divisor.NormFinite
import ABC3.Found.Divisor.GlobalUnits
import ABC3.Found.Divisor.SchemeOrdUnit
import ABC3.Found.FrdI.Sec6GaloisCat
import Mathlib.AlgebraicGeometry.Normalization
import Mathlib.AlgebraicGeometry.Morphisms.Proper

/-!
# 因子論の第 3 層 —— `V[L]`(相対正規化)と `D_L`(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> If we write V [L] for the normalization of V in L [so V [L] is also

原文 (FrdI p.109):
> for the set of prime divisors of V [L] that map into [possibly subvarieties of

## ★★この層の出口は `CartierDatum` である

底の圏 `𝒟 = B(G)⁰` は**閉じている**(`Found/FrdI/Sec6GaloisCat.lean` の
`(FinSub K K̄)ᵒᵖ` —— 連結・totally epimorphic・of FSM-type)。
単系論も**閉じている**(`Found/Divisor/Cartier*.lean`)。
組み立ても**閉じている**(`Found/Divisor/CartierFrobenioid.lean` の
`cartierFrobenioid` が `Theorem 5.2, (ii)` を当てる)。

★★**残っているのはこの層だけ**である —— `Found/Divisor/CartierFrobenioid.lean` の
`CartierDatum`(9 フィールド)を幾何から作ること。

## ★在庫測定(2026-08-20)

| 要るもの | mathlib | 判定 |
|---|---|---|
| 相対正規化 `f.normalization` | `AlgebraicGeometry/Normalization.lean` | ★★**ある**(`toNormalization` / `fromNormalization` つき) |
| `IsIntegralHom f.fromNormalization` | 同上(instance) | ★**ある** |
| proper 射 | `AlgebraicGeometry/Morphisms/Proper.lean` | ★**ある** |
| geometrically integral | `AlgebraicGeometry/Geometrically/Integral.lean` | ★**ある** |
| 「正規スキーム」の述語 | grep 0 件 | ★**無い**(`SchemeWeil.lean` で置いた) |
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta

universe u

/-! ## ★1. `V[L]` —— `L` における `V` の正規化(鎖 `normalize` の `normalization-in-ext2`) -/

variable {V : Scheme.{u}} {SpecL : Scheme.{u}}

/-- ★★**`V[L]`** —— `L` における `V` の正規化。

★mathlib の相対正規化 `f.normalization` がそのまま使える。
`f : Spec L ⟶ V` は函数体の有限拡大 `L/K` が定める射である。 -/
noncomputable abbrev normalizationIn (f : SpecL ⟶ V) [QuasiCompact f] [QuasiSeparated f] : Scheme.{u} := f.normalization

/-- ★★**`V[L]` も正規である**。

## ★★★★測定の訂正(2026-08-25)—— `V` の正規性は**要らなかった**

旧版は `_hV : IsNormalScheme V` を仮定に置いていたが、
**正規化は `V` が正規でなくても正規**である(整閉包は整閉)。
★仮定を落とした。

★★★★★**閉じた**(2026-08-25)—— 中身は `Found/Divisor/NormNormal.lean` にある。 -/
theorem isNormalScheme_normalizationIn {V : Scheme.{u}} [IsIntegral V] {Kbar : Type u} [Field Kbar]
    [Algebra V.functionField Kbar] (L : ABC3.Found.FrdI.FinSub V.functionField Kbar) :
    IsNormalScheme (normalizationIn (ABC3.Found.Divisor.specToV V L)) :=
  ABC3.Found.Divisor.normObj_isNormalScheme V L

/-- ★★**`V[L]` も proper である**(鎖 `normalize` の `normalization-proper-normal`)。

原文 (FrdI p.109):
> If we write V [L] for the normalization of V in L [so V [L] is also

★正規化射は整(mathlib の instance)で、有限拡大なら有限、したがって proper。
proper の合成は proper。

## ★★★★逸脱の記録(2026-08-25)—— 仮定が足りていなかった

旧版の仮定は `_hfin : LocallyOfFiniteType f`(`f : Spec L ⟶ V` の側)だったが、
**それでは `f.fromNormalization` の有限型性は出ない**。要るのは

* `Γ(V,U)` が**整閉かつ Noether**(原文の「proper normal variety」)
* `L/K` が**分離的**(原文の「necessarily separable」)

の 3 つで、そこから `IsIntegralClosure.finite`(mathlib)が
「整閉包は加群として有限」を与える。★仮定を原文に合わせて直した。

★★★★★**閉じた**(2026-08-25)—— 中身は `Found/Divisor/NormFinite.lean` にある。 -/
theorem isProper_normalizationIn {S V : Scheme.{u}} [IsIntegral V] {Kbar : Type u} [Field Kbar]
    [Algebra V.functionField Kbar] (L : ABC3.Found.FrdI.FinSub V.functionField Kbar)
    (g : V ⟶ S) (hg : IsProper g)
    (hnoeth : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsNoetherianRing Γ(V, U))
    (hic : ∀ U : V.Opens, IsAffineOpen U → (U : Set V).Nonempty → IsIntegrallyClosed Γ(V, U))
    (hsep : Algebra.IsSeparable V.functionField L.toIF) :
    IsProper (ABC3.Found.Divisor.normDown V L ≫ g) :=
  (ABC3.Found.Divisor.isProper_normObj V L g hg hnoeth hic hsep).1

/-! ## ★2. `D_L`(鎖 `normalize` の `DL-set`)

原文 (FrdI p.109):
> for the set of prime divisors of V [L] that map into [possibly subvarieties of

★「写る」は**像が `D_K` の素因子に含まれる**こと。
★原文が括弧で「possibly subvarieties of codimension ≥ 1 of」と断るとおり、
**余次元は上がってよい**(全射である必要はない)。 -/

/-- ★★**`D_L`** —— `D_K` の素因子に写る `V[L]` の素因子。 -/
def DLOf {W : Scheme.{u}} (π : W ⟶ V) (DK : Set (PrimeDivisorPt V)) :
    Set (PrimeDivisorPt W) :=
  {q | ∃ p ∈ DK, π.base q.1 ∈ closure ({p.1} : Set V)}

/-! ## ★3. `K`-`Q`-Cartier(鎖 `normalize` の `KQC`)

原文 (FrdI p.109):
> then we shall say that DK is K-Q-Cartier. In the following, we shall assume that

★原文はこれを**仮定**に置くので、我々も仮定として書く。 -/

/-- ★★**`K`-`Q`-Cartier** —— どの `Spec L` についても `D_L` の素因子が `Q`-Cartier。 -/
def IsKQCartier {ι : Type u} (W : ι → Scheme.{u})
    [∀ i, IsIntegral (W i)] [∀ i, AlgebraicGeometry.IsNoetherian (W i)]
    (DL : ∀ i, Set (PrimeDivisorPt (W i))) : Prop :=
  ∀ i, ∀ q ∈ DL i, IsQCartierDiv (W i) (Finsupp.single q 1)

/-! ## ★4. 大域函数(鎖 `normalize` の `global-units`)

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV, -/

open ABC3.Found.Divisor in
/-- ★★★★**代数的 Hartogs のスキーム版** —— `global-units` に残っている **1 本**。

★環の側(`Found/Divisor/Hartogs.lean` の `mem_Rsub_of_forall_heightOnePrime`)は
**閉じている**。残るのは配線だけで、

* `Γ(X,U)` の**高さ 1 素イデアル** ↔ `X` の**余次元 1 の点**(`isCodimOnePt_spec_iff`)
* `ordPt ≥ 0` ↔ 茎に入る(`exists_mem_stalk_of_ordPt_nonneg`、実装済み)
* アフィン開 `U ≅ Spec Γ(X,U)` を通す

の 3 段である。★貼り合わせの側(`exists_globalSection`)は**閉じた**。 -/
theorem mem_range_germ_of_forall_ordPt_nonneg {X : Scheme.{u}} [IsIntegral X]
    [AlgebraicGeometry.IsLocallyNoetherian X] (hnorm : IsNormalScheme X)
    (_hnoeth : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsNoetherianRing Γ(X, U))
    (u : X.functionField)
    (_h : ∀ x : PrimeDivisorPt X, 0 ≤ ordPt X hnorm x u) :
    u ∈ Set.range (CategoryTheory.ConcreteCategory.hom (X.germToFunctionField ⊤)) := by
  sorry

open ABC3.Found.Divisor in
/-- ★★**`O^×(A) = O^▷(A) = k_L^×`** —— proper normal なら大域函数は定数。

★`k_L` は `L` の中での `k` の代数閉包(`k` の有限分離拡大)。

## ★★★★測定の訂正(2026-08-25)—— 主張が偽だった

旧版は `Function.Bijective (g.appTop)`(`Γ(S) → Γ(W)` が全単射)と書いていたが、
**これは偽**である —— 原文が言うのは `Γ(W,⊤) = k_L` であって `k` ではない。
`k_L` は `k` の**有限拡大**で、一般に `k` とは違う。

★正しい形は「`ord` がすべての余次元 1 の点で `0` の元は `Γ(X,⊤)` の**単元**」であり、
そのうえで `Γ(X,⊤)` が `k` 上有限次の**体**(＝ `k_L`)である
(`Found/Divisor/GlobalUnits.lean` の `globalSections_isField` / `globalSections_finite`)。 -/
theorem globalSections_eq_constants {X : Scheme.{u}} [IsIntegral X]
    [AlgebraicGeometry.IsLocallyNoetherian X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    (hnorm : IsNormalScheme X)
    (hnoeth : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsNoetherianRing Γ(X, U))
    (u : X.functionField) (hu : u ≠ 0)
    (h : ∀ x : PrimeDivisorPt X, ordPt X hnorm x u = 0) :
    ∃ t : Γ(X, ⊤), IsUnit t ∧
      (CategoryTheory.ConcreteCategory.hom (X.germToFunctionField ⊤)) t = u :=
  exists_unit_globalSection g hg hu
    (mem_range_germ_of_forall_ordPt_nonneg hnorm hnoeth u (fun x => (h x).ge))

/-! ## ★5. 出口 —— `CartierDatum` を作る(鎖 `normalize` の `cartier-datum-geom`) -/

open ABC3.Found.FrdI in
/-- ★★★**幾何のデータから `CartierDatum` を作る**。

★★ここまで来ると `Found/Divisor/CartierFrobenioid.lean` の
`cartierFrobenioid`(`Theorem 5.2, (ii)` を当てたもの)が**そのまま**動き、
`Example 6.1` の model Frobenioid `C_{V,K̄,D_K}` が出る。

★入力は
* `VL L = V[L]`(相対正規化)
* `DL L ⊆ D_{V[L]}`(`D_K` の上にある素因子)
* `L ⟶ M` に対する `V[M] ⟶ V[L]`(底変換)
* `K`-`Q`-Cartier 性

★9 フィールドのうち `pull` / `pull_id` / `pull_comp` / `pull_mem` / `pull_nonneg` /
`pull_inj` は**すべて `Cartier.lean` の `pullbackCartier` の性質**である。 -/
theorem exists_cartierDatum_of_geometry
    (K Kbar : Type u) [Field K] [Field Kbar] [Algebra K Kbar]
    (VL : FinSub K Kbar → Scheme.{u})
    [inst1 : ∀ L, IsIntegral (VL L)] [inst2 : ∀ L, AlgebraicGeometry.IsNoetherian (VL L)]
    (_hnorm : ∀ L, IsNormalScheme (VL L))
    (DL : ∀ L, Set (PrimeDivisorPt (VL L)))
    (_base : ∀ {L M : FinSub K Kbar}, (L ⟶ M) → (VL M ⟶ VL L))
    (_hKQC : IsKQCartier VL DL) :
    Nonempty (CartierDatum.{u, u, u} ((FinSub K Kbar)ᵒᵖ)) := by
  sorry

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def normalizationIn.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — V[L](L における V の正規化)",
    sectionId := "frdi-example-6-1" }

def isNormalScheme_normalizationIn.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — V[L] は正規",
    sectionId := "frdi-example-6-1" }

def isNormalScheme_normalizationIn.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Found 側の本体(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.normObj_isNormalScheme") 109,
    .citation "[mathlib]" "Scheme.Hom.normalization(相対正規化)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalization") 109,
    .derivation "正規化は各アフィンで整閉包を取るので、茎は整閉" 109,
    .implicitStep "★原文は角括弧で「so V [L] is also a proper normal variety」と述べるだけ" 109 ]

def isProper_normalizationIn.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — V[L] は proper",
    sectionId := "frdi-example-6-1" }

def isProper_normalizationIn.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Found 側の本体(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isProper_normObj") 109,
    .citation "[mathlib]" "IsIntegralHom Scheme.Hom.fromNormalization(正規化射は整)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.fromNormalization") 109,
    .citation "[mathlib]" "IsIntegralClosure.finite(整閉包は加群として有限)"
      (.inMathlib "IsIntegralClosure.finite") 109,
    .citation "[mathlib]" "IsProper(proper 射と合成の安定性)"
      (.inMathlib "AlgebraicGeometry.IsProper") 109,
    .derivation "整閉包が加群として有限 ⟹ 正規化射は有限、有限射は proper、proper の合成は proper" 109,
    .implicitStep "★原文は角括弧で畳む" 109 ]

def DLOf.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — D_L(D_K の上にある素因子)",
    sectionId := "frdi-example-6-1" }

def IsKQCartier.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — K-Q-Cartier",
    sectionId := "frdi-example-6-1" }

def globalSections_eq_constants.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Example 6.1 — O^×(A) = O^▷(A) = k_L^×",
    sectionId := "frdi-example-6-1" }

def globalSections_eq_constants.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_unit_globalSection(Found、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_unit_globalSection") 110,
    .citation "[ABC3]" "mem_range_germ_of_forall_ordPt_nonneg(代数的 Hartogs のスキーム版、残り 1 本)"
      (.inProject "ABC3" "ABC3.Skeleton.Divisor.mem_range_germ_of_forall_ordPt_nonneg") 110,
    .citation "[mathlib]" "isField_of_universallyClosed(proper なら Γ(X,⊤) は体)"
      (.inMathlib "AlgebraicGeometry.isField_of_universallyClosed") 110,
    .citation "[mathlib]" "finite_appTop_of_universallyClosed(Γ(X,⊤) は k 上有限)"
      (.inMathlib "AlgebraicGeometry.finite_appTop_of_universallyClosed") 110,
    .derivation "proper + 整 ⟹ 大域切断は基礎体上有限次(したがって体、k の代数閉包)" 110,
    .implicitStep
      "★原文は「[since V [L] is a proper normal variety]」の角括弧ひとつで畳む" 110 ]

def mem_range_germ_of_forall_ordPt_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Example 6.1 — 代数的 Hartogs のスキーム版",
    sectionId := "frdi-example-6-1" }

/-- ★★★**`global-units` に残っている 1 本**。環の側は閉じているので、残るのは配線。 -/
def mem_range_germ_of_forall_ordPt_nonneg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "代数的 Hartogs(環の側、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.Hartogs.mem_Rsub_of_forall_heightOnePrime") 110,
    .citation "[ABC3]" "exists_globalSection(貼り合わせ、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_globalSection") 110,
    .citation "[ABC3]" "exists_mem_stalk_of_ordPt_nonneg(ord ≥ 0 なら茎に入る)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_mem_stalk_of_ordPt_nonneg") 110,
    .citation "[ABC3]" "isCodimOnePt_spec_iff(Spec R の余次元 1 の点は高さ 1 の素イデアル)"
      (.inProject "ABC3" "ABC3.Found.Divisor.isCodimOnePt_spec_iff") 110,
    .derivation
      "アフィン開 U ≅ Spec Γ(X,U) を通し、高さ 1 素イデアルごとに ord ≥ 0 を使って環の Hartogs を当て、貼り合わせる" 110 ]

def exists_cartierDatum_of_geometry.src : Source :=
  { paper := "FrdI", pdfPage := 109, item := "Example 6.1 — Φ・B が 𝒟 上の単系を定める",
    sectionId := "frdi-example-6-1" }

/-- ★★9 フィールドの内訳を書き下す。**引き戻し 1 本**に集約されている。 -/
def exists_cartierDatum_of_geometry.needs : List ProofObligation :=
  [ .citation "[ABC3]" "CartierDatum(9 フィールドの構造)"
      (.inProject "ABC3" "ABC3.Found.FrdI.CartierDatum") 109,
    .citation "[ABC3]" "pullbackCartier / pullbackCartier_add / isCartierDiv_pullbackCartier / pullbackCartier_nonneg"
      (.inProject "ABC3" "ABC3.Skeleton.Divisor.pullbackCartier") 109,
    .citation "[ABC3]" "isQCartierSubgroup_of_forall_isQCartier(qc フィールド)"
      (.inProject "ABC3" "ABC3.Skeleton.Divisor.isQCartierSubgroup_of_forall_isQCartier") 109,
    .citation "[ABC3]" "finSubOp_totallyEpimorphic / finSubOp_isOfFSMType(底の圏 B(G)⁰)"
      (.inProject "ABC3" "ABC3.Found.FrdI.finSubOp_isOfFSMType") 109,
    .derivation
      "pull_id / pull_comp は引き戻しの関手性、pull_inj は V[M] → V[L] が支配的であることから" 109,
    .implicitStep
      "★原文は「the assignments L → Φ(L), L → B(L) determine, respectively, a perf-factorial divisorial monoid Φ on D」と一文で述べる" 109 ]

end ABC3.Skeleton.Divisor
