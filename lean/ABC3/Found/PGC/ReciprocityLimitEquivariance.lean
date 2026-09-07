import ABC3.Found.PGC.SemilinearRestriction

/-!
# 極限への持ち上げと、`E` を `ρ` 由来のものに取り替える段

直前の節点(`SemilinearRestriction.lean`)は「cross-point instance bridging」を
正面突破して **(ii) の有限段** `ρ_{f^σ,n}(Φ τ Φ^{-1}) = σ(ρ_{f,n}(τ))`
(`reciprocityMap_semilinear_conj`)を出し、残りを 4 つ名指しした:

> 1. 極限への持ち上げ(`reciprocityMap` の**点非依存性**が要る。木には無い)。
> 2. `Gal(L^{ab}/L) ≅ 𝒪_L^× × Ẑ` の `E` を `ρ` 由来のものに取り替える段(★次の壁)。
> 3. (iii) 素元非依存(2 の内側で要る)。
> 4. `Γ_F` への代入(機械的)。

本ファイルは **1 と 4 を完全に埋め**、**2 を(3 を除いて)全部埋めた**。
その結果 `ArtinUnitEquivariance`(壁)は **1 本の名前の付いた穴**
(`ReciprocityDatumIndependenceOnTorsion` ＝ 原典 Corollary 4.9 を `p` 冪捩れの
上に制限したもの)に**還元された**。

## ★★何がこのファイルで出たか(測定結果)

| 段 | 宣言 | 状態 |
|---|---|---|
| 核 A | `act_cocycle_indep` | ★証明した —— 単純推移的な作用の上のコサイクルは基点に依らない(純代数) |
| 核 B | `sub_mem_span_pow_map` | ★証明した —— 環準同型は合同を合同へ(純環論) |
| 1-(a) | `reciprocityMap_point_indep` | ★★★★証明した —— `ρ_{f,n}` は評価点に依らない(木に無かった) |
| 1 | `reciprocityUnits_semilinear_conj` | ★★★★★証明した —— (ii) の極限段 |
| 4 | `reciprocityUnits_absGalConjCME` | ★★証明した —— `Γ_F` への代入 |
| 2-(a) | `abelianGalEquivUnitsZHat_restrictNormalHom_fst` | ★★★証明した —— `E` の単数成分は `ρ` そのもの |
| 2-(b) | `abelianGalConj_restrictNormalHom` | ★★証明した —— `Ψ_g` と `Gal(L^{ab}/L)` の共役の両立 |
| 2 | `artinUnitEquivariance_of_reciprocityDatumIndependenceOnTorsion` | ★★★★★証明した(穴 1 本を仮定して) |
| 3 | `ReciprocityDatumIndependenceOnTorsion` | ★★残った(仮定として名前を付けた) |
| 壁 | `ArtinUnitEquivariance` | ★★**仮定つきで出た。無条件ではまだ出ていない。** |

★★**`ArtinUnitEquivariance` は無条件には出ていない。**出ていないものを出たと書かない。

## ★★★1 はどう埋まったか —— 「点非依存性」は純代数だった

申し送りは「`reciprocityMap` の点非依存性が木に無い(`reciprocityMap_congr` は
`x = x'` の場合だけ)。本ファイルの §10(点をまたぐ同変性)+ 作用の乗法性で
出せる見込み。★測っていない」と書いていた。★**測ったら、そのとおりに出た。**

抽象核はこうである(★分岐・付値・Galois の語彙が 1 語も出てこない):

> 可換群 `U` が集合 `X` に作用し、`F : X → X` がその作用と可換で、
> 作用が基点 `x` で自由なら、「`F x = u · x` を満たす `u`」は基点に依らない。

証明は 4 行:`F x' = u'·(w·x) = (u'w)·x` と `F x' = w·(F x) = (wu)·x` を比べて
自由性で `u'w = wu`、可換性で `u = u'`。
★具体層では `U := 𝒪^×/(1+π^n)`、`X := ψ_n の根`、`F := σ` を代入する。
* 乗法性 `(ab)·y = a·(b·y)` は在庫(`lubinTateActionAtTorsionPoint_mul_eq_of_action`)。
* 同変性 `σ(a·y) = a·(σ y)` は**直前の波の副産物**(`algEquiv_lubinTateActionAtTorsionPoint_cross`)。
* 自由性は在庫(`unitActionQuotientLift_injective`)。
* 推移性(`x' = w·x`)は在庫(`unitActionQuotientBijOn_bijective`)。

★★**したがって「点非依存性」に新しい数学は 1 つも要らなかった** ——
直前の波が §10 で cross-point を突破した時点で、材料は全部そろっていた。

## ★★極限への持ち上げの中身

`reciprocityUnits = (𝒪^× ≅ lim (𝒪/π^n)^×)^{-1} ∘ ρ_limit` なので、
単数の等式は**各 `n` での合同**に落ちる:

* `n = 0`:`(𝒪/π^0)^×` は自明群。
* `n = k+1`:族の第 `k+1` 成分は `ρ_{f,k+1}` を生成元 `x_k`(`psiGenSeq`)の上で
  評価したもの。★**捩れた塔の生成元 `y_k` は `Φ x_k` とは限らない**が、
  点非依存性で移せる。あとは有限段(`reciprocityMap_semilinear_conj`)と、
  「合同は環準同型で保たれる」(核 B)だけ。

★★**申し送りの「生成元列の付け替え」は、点非依存性 1 本に吸収された。**

## ★★★2 はどこまで行けたか(★次の壁と名指しされていた段)

申し送りは
> `nonempty_abelianGalContinuousEquivUnitsZHat` が与えるのは同型の存在だけで、
> それが `ρ` から来ていることを言っていない。
と書いていた。★**測ったら、`ρ` 由来であることは木の在庫から 3 行で出た。**

`abelianGalEquivProd_restrictNormalHom`(`ArtinMap.lean`、Λ8)が
`Φ(τ|_M) = (ρ(τ|_{K_π}), τ|_{K^ur})` を与え、
`lubinTateClosureGalEquivUnits_restrictNormalHom`(`LubinTateClosureTopology.lean`)が
`ρ(τ|_{K_π}) = reciprocityUnits τ` を与える。
★**つまり「存在だけ」だったのは `Nonempty` に包んだ形の方だけで、
包む前の `abelianGalEquivUnitsZHat` は最初から `ρ` 由来だった。**

残ったのは **3(素元非依存)だけ**である。しかも要るのは
★**`p` 冪捩れの上だけ**(`Ẑ` 成分は捩れが無いので、捩れ元は自動的に慣性に入る)。
これを `ReciprocityDatumIndependenceOnTorsion` と名付けて仮定に立てた。

## ★★残った穴が「本当に必要な分だけ」であることの確認

`ReciprocityDatumIndependenceOnTorsion` は
「**同じ `K` の上の 2 つの Lubin-Tate データ `(π,f)`・`(π',f')` について、
`Gal(K^{ab}/K)` への制限が `p^m` 捩れである `θ` の上では `ρ_f(θ) = ρ_{f'}(θ)`**」
という主張である。これは原典 Corollary 4.9(「`ρ_f` は `f` に依らない」)を
**捩れの上に制限したもの**であり、古典的に正しい:
`Gal(K^{ab}/K) ≅ 𝒪_K^× × Ẑ` の `Ẑ` は捩れを持たないので、捩れ元は慣性部分に入り、
慣性の上では Artin 写像は素元の選び方に依らない。

★★**捩れの条件を落とすと偽になる**(直前の波の測定どおり `ρ_f ≠ ρ_{f^σ}`)。
例:`θ := Art_{π'}(π')` は `ρ_{f'}(θ) = 1` だが `ρ_f(θ) = π'/π ≠ 1`。
★**だから条件付きの形にした。通すために偽の仮定は置いていない。**

## ★★段取りとの差分

* ★★**`LubinTateEndoTwisted` は要らなかった**(本体の見立ては 5 波連続で外れ)。
  ★**`powerSeries_uniqueness` も出てこない。**
* ★★**素元非依存性(iii)は「使わなかった」のではなく「ここで初めて要った」。**
  ★ただし要ったのは**捩れの上だけ**で、しかも**仮定として名前を付けられる形**だった。
  ★申し送りの見立て(「慣性部分＝`p^n` 捩れの上だけでよい見込み」)は当たっていた。
* ★★**`reciprocityMap` の点非依存性は「新しい数学」ではなかった** ——
  直前の波の §10 の副産物を抽象核に流し込むだけで出た(一発)。
* ★★**`E` が `ρ` 由来であることは既に木に在った。**「壁」の半分は在庫だった。
* ★**#59 には 1 度も当たっていない。** 中間体は `K⟮x⟯` の 1 層だけで、
  `restrictScalars`(体の)も書いていない。`restrictNormalHom` は**書いたが 1 層**
  (`K̄ → K^{ab}`)で、中間体の中の中間体は作っていない(#59 の定型 (b))。
* ★**#211 に当たった**(予告どおり)。`(fixedFieldAut S g).restrictScalars ℚ_[p]` を
  `integerRingEquiv` の裸の引数に置くと `CoeT` の探索が timeout する。
  ★**在庫の `fixedFieldIntegerAut`(型注釈つきの `def`)に差し替えて回避した。**

## ★★退化の自己検査

* ★**`E` の存在だけでは足りない**ので、`abelianGalEquivUnitsZHat`(構成的な方)を使い、
  その単数成分が `ρ` であることを**定理として**述べた。
* ★**`ArtinUnitEquivariance` は `p^n` 捩れの上でしか要求していない**ので、
  穴もそこだけに限った(`Ẑ` 成分は条件に現れない)。
* ★**`Φ x` は `K⟮x⟯` に留まらない。** 本ファイルもどこでも留まると仮定していない。
* ★**易しい半分(`g ∈ S`)は壊していない** —— `ArtinUnitEquivariance` の定義自体は
  1 文字も触っていない。
* ★`ℕ∞` の切り詰め引き算・除算は書いていない(#102)。

## ★★残っているものを名指しで

1. ★★★**`ReciprocityDatumIndependenceOnTorsion`** ——
   原典 Corollary 4.9 の本体(`θ ∈ Θ^{K̂}_{π,π'}` を Dwork/Lubin-Tate で作り、
   `reciprocity_eq_of_intertwiner`(`LubinTateReciprocityIndependence.lean`)に代入する)。
   ★**捩れの上に制限してあるので、`Ẑ` 成分の議論は要らない。**
   ★在庫:`LubinTateUniformizerIndependence.lean` / `DworkMultiplicative.lean` /
   `reciprocity_eq_of_thetaSet`。
2. これが埋まれば `ArtinUnitEquivariance` は**無条件で出る**
   (本ファイルの `artinUnitEquivariance_of_reciprocityDatumIndependenceOnTorsion`)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

open scoped NormedField Valued Classical

/-! ## 1. 抽象核 A —— 単純推移的な作用の上のコサイクルは基点に依らない

★この節に分岐・付値・Galois・Lubin-Tate の語彙は 1 つも出てこない。
可換群 `U` の集合 `X` への作用と、それと可換な写像 `F : X → X` の話である。

★`MulAction` インスタンスを**要求しない**形にしてある(`act` を裸の関数で受ける)。
具体層の作用は「点ごとに証明項を渡して定義する」形なので、インスタンスに
束ねると却って高くつく(`LubinTateReciprocityIndependence.lean` の §1b と同じ判断)。 -/

section AbstractCoreA

variable {U X : Type*} [CommGroup U]

/-- ★★★★**抽象核 A** —— `F x = u · x` を満たす `u` は基点 `x` に依らない。

仮定は 4 つだけ:
* `act` は結合的(`a · (b · y) = (ab) · y`)、
* `F` は作用と可換(`F (a · y) = a · (F y)`)、
* 基点 `x` から `x'` へ移す `w` がある(推移性)、
* 基点 `x` で作用は自由。

★★証明は 4 行:`F x' = u'·(w·x) = (u'w)·x` と `F x' = w·(F x) = (wu)·x` を
比べて自由性で `u'w = wu`、`U` の可換性で `u = u'`。

★これが `reciprocityMap` の**点非依存性**の中身の全部である。 -/
theorem act_cocycle_indep (act : U → X → X) (F : X → X)
    (hmul : ∀ (a b : U) (y : X), act a (act b y) = act (a * b) y)
    (hF : ∀ (a : U) (y : X), F (act a y) = act a (F y))
    {x x' : X} {w u u' : U}
    (hw : act w x = x')
    (hfree : ∀ a b : U, act a x = act b x → a = b)
    (hu : act u x = F x) (hu' : act u' x' = F x') : u = u' := by
  have h1 : act (u' * w) x = F x' := by rw [← hmul, hw, hu']
  have h2 : act (w * u) x = F x' := by rw [← hmul, hu, ← hF, hw]
  have h3 : u' * w = w * u := hfree _ _ (h1.trans h2.symm)
  rw [mul_comm w u] at h3
  exact (mul_right_cancel h3).symm

end AbstractCoreA

variable {p : ℕ} [Fact p.Prime]

/-! ## 2. 具体層 1 —— `ψ_n` の根の集合の上の作用 -/

open scoped Classical in
/-- `ψ_n` の根は `Λ_n` の元 —— 在庫の `Λ_n \ Λ_{n-1} = (ψ_n の根)` から。 -/
theorem mem_iteratedLubinTateTorsionPoints_of_mem_psi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) {y : K.closure}
    (hy : y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    y ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn] at hy
  exact (Finset.mem_sdiff.mp hy).1

section PointIndependence

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n)

open scoped Classical in
/-- `ψ_n` の根の集合の上の `𝒪^×/(1+π^n)` 作用。 -/
noncomputable def psiRootAct
    (U : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)
    (y : (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure)) :
    (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure) :=
  ⟨(↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n (y : K.closure)
        (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
        (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure)) U) :
      IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) : K.closure),
    unitActionQuotientLift_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn
      (y : K.closure) y.2
      (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
      (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure)) U⟩

open scoped Classical in
@[simp] theorem coe_psiRootAct
    (U : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)
    (y : (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure)) :
    ((psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn U y : K.closure))
      = (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n (y : K.closure)
          (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
          (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure)) U) :
        IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) : K.closure) := rfl

open scoped Classical in
/-- 作用の結合性。 -/
theorem psiRootAct_mul
    (a b : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)
    (y : (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure)) :
    psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn a
        (psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn b y)
      = psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn (a * b) y := by
  induction a using QuotientGroup.induction_on with
  | H u =>
  induction b using QuotientGroup.induction_on with
  | H v =>
  apply Subtype.ext
  exact (lubinTateActionAtTorsionPoint_mul_eq_of_action K hq hπmax hπne0 f hf0 hf1 hf n
    (y : K.closure)
    (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
    (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure))
    (v : 𝒪[K.carrier]) (u : 𝒪[K.carrier])
    (lubinTateActionAtTorsionPoint_mem K hq hπmax hπne0 f hf0 hf1 hf n (y : K.closure)
      (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
      (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure)) (v : 𝒪[K.carrier]))
    (IntermediateField.mem_adjoin_simple_self K.carrier _)).symm

open scoped Classical in
/-- Galois の作用を `ψ_n` の根の集合の上の写像として見たもの。 -/
noncomputable def psiRootGal (σ : K.closure ≃ₐ[K.carrier] K.closure)
    (y : (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure)) :
    (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure) :=
  ⟨σ (y : K.closure), algEquiv_mem_iteratedLubinTatePsiTorsionPoints_of_mem
    K hq hπmax hπne0 f hf0 hf1 hf n hn σ (y : K.closure) y.2⟩

open scoped Classical in
/-- `σ(a·y) = a·(σ y)`。 -/
theorem psiRootGal_act (σ : K.closure ≃ₐ[K.carrier] K.closure)
    (a : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)
    (y : (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure)) :
    psiRootGal K hq hπmax hπne0 f hf0 hf1 hf n hn σ
        (psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn a y)
      = psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn a
          (psiRootGal K hq hπmax hπne0 f hf0 hf1 hf n hn σ y) := by
  induction a using QuotientGroup.induction_on with
  | H u =>
  apply Subtype.ext
  exact algEquiv_lubinTateActionAtTorsionPoint_cross K hq hπmax hπne0 f hf0 hf1 hf n
    (y : K.closure)
    (mem_iteratedLubinTateTorsionPoints_of_mem_psi K hq hπmax hπne0 f hf0 hf1 hf n hn y.2)
    (IntermediateField.mem_adjoin_simple_self K.carrier (y : K.closure)) σ (u : 𝒪[K.carrier])

/-! ## 3. ★★★★払い出し A —— `ρ_{f,n}` の点非依存性 -/

def reciprocityMap_point_indep.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 8, item := "Proposition 4.7", sectionId := "prop-4-7" }

open scoped Classical in
/-- ★★★★**`ρ_{f,n}` は評価点の取り方に依らない**。

原文 (Yoshida08 p.8):
> The L^m_f is Galois over K, and the following map is bijective for any α ∈ µ^×_f,m:

★原典は「任意の `α`(原始的な `π^m` 捩れ点)について同じ全単射が得られる」と
述べており、その「任意の `α` について」の部分がここで形式化された。
木の `reciprocityMap_congr` は `x = x'` の場合(証明項だけが違う場合)しか無く、
★**本当に異なる 2 点をまたぐ主張は無かった。**

★★証明は抽象核 A に
* 推移性 = `unitActionQuotientBijOn_bijective`(在庫)、
* 自由性 = `unitActionQuotientBijOn_injective`(在庫)、
* 乗法性 = `psiRootAct_mul`(在庫 `lubinTateActionAtTorsionPoint_mul_eq_of_action` の言い換え)、
* 同変性 = `psiRootGal_act`(直前の波の `algEquiv_lubinTateActionAtTorsionPoint_cross`)

を代入するだけである。★**新しい数学は 1 つも要らなかった。** -/
theorem reciprocityMap_point_indep (x x' : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hx'ψ : x' ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hx'n : x' ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem' : x' ∈ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
      = reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x' hx'ψ hx'n hmem' σ := by
  obtain ⟨W, hW⟩ := (unitActionQuotientBijOn_bijective K hq hπmax hπne0 f hf0 hf1 hf n hn
    x hxψ hxn hmem).2 ⟨x', hx'ψ⟩
  refine act_cocycle_indep (psiRootAct K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (psiRootGal K hq hπmax hπne0 f hf0 hf1 hf n hn σ)
    (psiRootAct_mul K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (psiRootGal_act K hq hπmax hπne0 f hf0 hf1 hf n hn σ)
    (x := ⟨x, hxψ⟩) (x' := ⟨x', hx'ψ⟩) (w := W) hW
    (fun a b h => unitActionQuotientBijOn_injective K hq hπmax hπne0 f hf0 hf1 hf n hn
      x hxψ hxn hmem h)
    (Subtype.ext (reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ))
    (Subtype.ext (reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x' hx'ψ hx'n hmem' σ))

end PointIndependence

/-! ## 4. 抽象核 B —— 合同は環準同型で保たれる

★この節も純粋な環論である(体・付値・Galois は出てこない)。 -/

section AbstractCoreB

/-- ★**抽象核 B** —— 環準同型は「`π^n` を法とする合同」を「`φπ^n` を法とする合同」に移す。 -/
theorem sub_mem_span_pow_map {A B : Type*} [CommRing A] [CommRing B] (φ : A →+* B)
    {π : A} {n : ℕ} {a b : A} (h : a - b ∈ Ideal.span ({π ^ n} : Set A)) :
    φ a - φ b ∈ Ideal.span ({(φ π) ^ n} : Set B) := by
  rw [← map_sub]
  rw [Ideal.mem_span_singleton] at h ⊢
  obtain ⟨c, hc⟩ := h
  exact ⟨φ c, by rw [hc, map_mul, map_pow]⟩

end AbstractCoreB

/-! ## 5. ★★★★★払い出し B —— 極限への持ち上げ(申し送りの 1) -/

section LimitLift

/-- `(𝒪_K/π^0)^×` は自明群(`π^0 = 1` なので剰余環が自明)。 -/
theorem subsingleton_units_quotient_pow_zero (K : PAdicLocalField p) (π : 𝒪[K.carrier]) :
    Subsingleton (𝒪[K.carrier] ⧸ Ideal.span ({π ^ 0} : Set 𝒪[K.carrier]))ˣ := by
  have h0 : Ideal.span ({π ^ 0} : Set 𝒪[K.carrier]) = ⊤ := by
    rw [pow_zero]; exact Ideal.span_singleton_one
  rw [h0]
  constructor
  intro a b
  apply Units.ext
  exact Subsingleton.elim _ _

/-- `reciprocityUnits` の第 `m` 成分は `reciprocityMapLimitFamily` の第 `m` 成分。 -/
theorem unitReductionQuotientMap_reciprocityUnits (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (τ : K.closure ≃ₐ[K.carrier] K.closure) (m : ℕ) :
    unitReductionQuotientMap K π m (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ)
      = reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf τ m :=
  congrFun (congrArg Subtype.val
    (MulEquiv.apply_symm_apply (unitsEquivCompatibleUnits K hπmax)
      (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf τ))) m

def reciprocityUnits_semilinear_conj.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

set_option maxHeartbeats 1000000 in
open scoped Classical in
/-- ★★★★★★**申し送りの 1(極限への持ち上げ)そのもの**:

  `ρ_{f^σ}(Φ τ Φ^{-1}) = σ(ρ_f(τ))` (単数群 `𝒪_K^×` の中の等式)。

原文 (Yoshida08 p.9):
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

★★**原典と同じ主張ではない**(逸脱として記録する)。原典 Corollary 4.9 は
「`ρ_f` は `f` に依らない」(同じ `K` の上の 2 つの `f` を比べる)であるのに対し、
ここでは `σ`(体の半線型自己同型)で捩った塔との比較で、係数側にも `σ` が現れる。
★役割(「極限の `ρ` が `f` の取り替えの下でどう動くか」の段)が同じなので
`.src` を Corollary 4.9 に付けた。

★★**証明の骨**:`reciprocityUnits` は `𝒪_K^× ≅ lim_n (𝒪_K/π^n)^×` の逆写像を
`ρ_limit` に合成したものなので、単数の等式は各 `n` での合同に落ちる。
* `n = 0`:自明群(`subsingleton_units_quotient_pow_zero`)。
* `n = k+1`:族の第 `k+1` 成分は `ρ_{f,k+1}` を `psiGenSeq` の生成元 `x_k` の上で
  評価したもの。★捩れた塔の生成元 `y_k` は `Φ x_k` とは限らないが、
  **点非依存性**(`reciprocityMap_point_indep`)で移せる。あとは有限段
  (`reciprocityMap_semilinear_conj`)と抽象核 B(合同の輸送)だけ。

★★**申し送りの「3. 生成元列の付け替え」は、点非依存性 1 本に吸収された。** -/
theorem reciprocityUnits_semilinear_conj (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    (σ : K.carrier ≃ₐ[ℚ_[p]] K.carrier) (Φ : K.closure ≃+* K.closure)
    (hΦ : ∀ a : K.carrier, Φ (algebraMap K.carrier K.closure a)
      = algebraMap K.carrier K.closure (σ a))
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (τ : K.closure ≃ₐ[K.carrier] K.closure) :
    ((reciprocityUnits K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf)
        (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ) : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])
      = integerRingEquiv σ ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ :
          (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) := by
  have hmain : reciprocityUnits K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf)
        (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ)
      = Units.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier])
          (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ) := by
    apply (unitsEquivCompatibleUnits K
      (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)).injective
    apply Subtype.ext
    funext m
    show unitReductionQuotientMap K ((integerRingEquiv σ) π) m _
      = unitReductionQuotientMap K ((integerRingEquiv σ) π) m _
    cases m with
    | zero => exact @Subsingleton.elim _ (subsingleton_units_quotient_pow_zero K _) _ _
    | succ k =>
      obtain ⟨u, hu⟩ := QuotientGroup.mk_surjective
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hmem τ)
      have hsub : (u : 𝒪[K.carrier])
          - ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ : (𝒪[K.carrier])ˣ) :
              𝒪[K.carrier]) ∈ Ideal.span ({π ^ (k + 1)} : Set 𝒪[K.carrier]) := by
        rw [← Ideal.Quotient.mk_eq_mk_iff_sub_mem]
        have h1 : unitReductionQuotientMap K π (k + 1) u
            = unitReductionQuotientMap K π (k + 1)
                (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ) := by
          rw [unitReductionQuotientMap_reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ (k + 1),
            ← principalUnitsQuotientEquiv_apply_mk K hπmax (k + 1) (by omega) u, hu]
          rfl
        exact congrArg Units.val h1
      have hpt := reciprocityMap_point_indep K hq
        (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
        (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
        (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (integerRingEquiv σ) hf) (k + 1) (by omega)
        (psiGenSeq K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
          (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
          (map_residue_map_twist (integerRingEquiv σ) hf) k).pt
        (Φ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt)
        (psiGenSeq K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
          (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
          (map_residue_map_twist (integerRingEquiv σ) hf) k).hψ
        (psiGenSeq K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
          (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
          (map_residue_map_twist (integerRingEquiv σ) hf) k).hn
        (psiGenSeq K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
          (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
          (map_residue_map_twist (integerRingEquiv σ) hf) k).hmem
        (map_mem_iteratedLubinTatePsiTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf (k + 1) (by omega) (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ)
        (map_mem_iteratedLubinTateTorsionPoints K hq (integerRingEquiv σ) hπmax hπne0
          f hf0 hf1 hf (k + 1) (Φ : K.closure →+* K.closure)
          (semilinear_integerRingEquiv K σ (Φ : K.closure →+* K.closure) hΦ)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn)
        (IntermediateField.mem_adjoin_simple_self K.carrier
          (Φ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt))
        (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ)
      have hconj := reciprocityMap_semilinear_conj K hq σ Φ hΦ hπmax hπne0 f hf0 hf1 hf
        (k + 1) (by omega) (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hmem τ u
        (Units.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) u) hu rfl
      have hstep : unitReductionQuotientMap K ((integerRingEquiv σ) π) (k + 1)
            (reciprocityUnits K hq (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
              (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
              (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
                𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
              (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
              (map_residue_map_twist (integerRingEquiv σ) hf)
              (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ))
          = unitReductionQuotientMap K ((integerRingEquiv σ) π) (k + 1)
              (Units.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
                𝒪[K.carrier] →+* 𝒪[K.carrier]) u) := by
        refine Eq.trans (unitReductionQuotientMap_reciprocityUnits K hq
          (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax)
          (fun h => hπne0 ((integerRingEquiv σ).injective (by rw [h, map_zero])))
          (PowerSeries.map ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
            𝒪[K.carrier] →+* 𝒪[K.carrier]) f)
          (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
          (map_residue_map_twist (integerRingEquiv σ) hf)
          (conjSemilinearAlgEquiv Φ σ.toRingEquiv hΦ τ) (k + 1)) ?_
        refine Eq.trans (congrArg (principalUnitsQuotientEquiv K
          (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax) (k + 1) (by omega))
          (hpt.trans hconj)) ?_
        exact principalUnitsQuotientEquiv_apply_mk K
          (maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax) (k + 1) (by omega) _
      refine hstep.trans ?_
      apply Units.ext
      show Ideal.Quotient.mk (Ideal.span ({((integerRingEquiv σ) π) ^ (k + 1)} :
          Set 𝒪[K.carrier])) (integerRingEquiv σ (u : 𝒪[K.carrier]))
        = Ideal.Quotient.mk (Ideal.span ({((integerRingEquiv σ) π) ^ (k + 1)} :
            Set 𝒪[K.carrier])) (integerRingEquiv σ
              ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ : (𝒪[K.carrier])ˣ) :
                𝒪[K.carrier]))
      rw [Ideal.Quotient.mk_eq_mk_iff_sub_mem]
      exact sub_mem_span_pow_map
        ((integerRingEquiv σ : 𝒪[K.carrier] ≃+* 𝒪[K.carrier]) :
          𝒪[K.carrier] →+* 𝒪[K.carrier]) hsub
  exact congrArg (fun w : (𝒪[K.carrier])ˣ => (w : 𝒪[K.carrier])) hmain

end LimitLift

/-! ## 6. ★★★払い出し C —— `Γ_F` の層への代入(申し送りの 4) -/

section GammaF

open scoped Classical in
/-- ★★★**`Γ_F` の層への代入**:`L := L_S` の上で

  `ρ_{f^{σ_g}}(Ψ_g τ) = σ_g(ρ_f(τ))`。

★直前の波の `absGalConjCME_eq_conjSemilinearAlgEquiv`(`Ψ_g` は抽象核 C の共役
そのもの、`rfl`)で書き換えて、上の極限段に代入するだけ(2 行)。

★★**#211 に当たった**:`(fixedFieldAut S g).restrictScalars ℚ_[p]` を
`integerRingEquiv` の裸の引数に置くと `CoeT` の探索が timeout する
(`K'` が未知のままなので `PAdicLocalField.carrier ?K'` への coe を探しに行く)。
★**在庫の `fixedFieldIntegerAut`(型注釈つきの `def`)に差し替えて回避した。** -/
theorem reciprocityUnits_absGalConjCME (F : PAdicLocalField p) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal))
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[(fixedFieldLocalField F S hopen).carrier])
      𝒪[(fixedFieldLocalField F S hopen).carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField
      𝒪[(fixedFieldLocalField F S hopen).carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[(fixedFieldLocalField F S hopen).carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField
      𝒪[(fixedFieldLocalField F S hopen).carrier]) = pp ^ ff) (g : F.absGal)
    {π : 𝒪[(fixedFieldLocalField F S hopen).carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[(fixedFieldLocalField F S hopen).carrier]
      = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[(fixedFieldLocalField F S hopen).carrier])
    (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue
      𝒪[(fixedFieldLocalField F S hopen).carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (τ : (fixedFieldLocalField F S hopen).absGal) :
    ((reciprocityUnits (fixedFieldLocalField F S hopen) hq
        (maximalIdeal_eq_span_map (fixedFieldIntegerAut F S hopen g) hπmax)
        (fun h => hπne0 ((fixedFieldIntegerAut F S hopen g).injective (by rw [h, map_zero])))
        (PowerSeries.map ((fixedFieldIntegerAut F S hopen g :
            𝒪[(fixedFieldLocalField F S hopen).carrier] ≃+*
              𝒪[(fixedFieldLocalField F S hopen).carrier]) :
          𝒪[(fixedFieldLocalField F S hopen).carrier] →+*
            𝒪[(fixedFieldLocalField F S hopen).carrier]) f)
        (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
        (map_residue_map_twist (fixedFieldIntegerAut F S hopen g) hf)
        (absGalConjCME F S hopen g τ) :
        (𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ) :
        𝒪[(fixedFieldLocalField F S hopen).carrier])
      = fixedFieldIntegerAut F S hopen g
          ((reciprocityUnits (fixedFieldLocalField F S hopen) hq hπmax hπne0 f hf0 hf1 hf τ :
            (𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ) :
            𝒪[(fixedFieldLocalField F S hopen).carrier]) := by
  rw [absGalConjCME_eq_conjSemilinearAlgEquiv]
  exact reciprocityUnits_semilinear_conj (fixedFieldLocalField F S hopen) hq
    ((fixedFieldAut S g).restrictScalars ℚ_[p]) (fixedFieldClosureAut F S hopen g)
    (fixedFieldClosureAut_semilinear F S hopen g) hπmax hπne0 f hf0 hf1 hf τ

end GammaF

/-! ## 7. ★★★★★申し送りの 2 —— `E` を `ρ` 由来のものに取り替える段

★申し送りは「`nonempty_abelianGalContinuousEquivUnitsZHat` が与えるのは同型の
存在だけで、それが `ρ` から来ていることを言っていない」を**次の壁**と名指しした。
★★**測ったら、`Nonempty` に包む前の `abelianGalEquivUnitsZHat` は最初から
`ρ` 由来で、そのことは木の在庫 2 本から出た。** -/

section RhoDerived

open scoped Classical in
/-- ★★★★**`E` の単数成分は `ρ` そのもの**:

  `(E (τ|_{K^{ab}})).1 = reciprocityUnits τ`。

★★これが申し送りの 2 の本体である。証明は在庫 2 本を継ぐだけ:
* `abelianGalEquivProd_restrictNormalHom`(`ArtinMap.lean`):
  `Φ(τ|_M) = (ρ(τ|_{K_π}), τ|_{K^{ur}})`、
* `lubinTateClosureGalEquivUnits_restrictNormalHom`(`LubinTateClosureTopology.lean`):
  `ρ(τ|_{K_π}) = reciprocityUnits τ`。

★`galEquivUnitsZHat` は第 2 成分だけを `Ẑ` に取り替えるので、第 1 成分は素通し。 -/
theorem abelianGalEquivUnitsZHat_restrictNormalHom_fst (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (τ : K.closure ≃ₐ[K.carrier] K.closure) :
    (abelianGalEquivUnitsZHat K hq hπmax hπne0 f hf0 hf1 hf
        (AlgEquiv.restrictNormalHom (abelianClosure K) τ)).1
      = reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf τ := by
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  rw [abelianGalEquivUnitsZHat, galEquivUnitsZHat_apply,
    abelianGalEquivProd_restrictNormalHom, lubinTateClosureGalEquivUnits_restrictNormalHom]

/-- ★★**`Ψ_g` の共役は `Gal(L^{ab}/L)` の共役の持ち上げ**:

  `abelianGalConj g (τ|_{L^{ab}}) = (Ψ_g τ)|_{L^{ab}}`。

★`abelianGalConj` は `Γ_L^{ab}(位相的)≅ Gal(L^{ab}/L)` を経由して定義されている
ので、`topAbelianizationEquivAbelianGal_mk`(値は制限写像)と
`topAbelianizationCME_mk`(`rfl`)の 2 本で `mk` の上へ降ろすだけ。 -/
theorem abelianGalConj_restrictNormalHom (F : PAdicLocalField p) (S : Subgroup F.absGal)
    [S.Normal] (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal)
    (τ : (fixedFieldLocalField F S hopen).absGal) :
    abelianGalConj F S hopen g
        (AlgEquiv.restrictNormalHom (abelianClosure (fixedFieldLocalField F S hopen)) τ)
      = AlgEquiv.restrictNormalHom (abelianClosure (fixedFieldLocalField F S hopen))
          (absGalConjCME F S hopen g τ) := by
  rw [abelianGalConj_apply, ← topAbelianizationEquivAbelianGal_mk, MulEquiv.symm_apply_apply,
    topAbelianizationCME_mk, topAbelianizationEquivAbelianGal_mk]

end RhoDerived

/-! ## 8. ★★★★★壁への還元 —— 残った穴は 1 本 -/

section ArtinReduction

def ReciprocityDatumIndependenceOnTorsion.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

/-- ★★★★★**残った穴** —— 原典 Corollary 4.9 を `p` 冪捩れの上に制限したもの。

原文 (Yoshida08 p.9):
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

主張:同じ `K` の上の 2 つの Lubin-Tate データ `(π,f)`・`(π',f')` について、
`Gal(K^{ab}/K)` への制限が `p^m` 捩れであるような `θ` の上では
`ρ_{f'}(θ) = ρ_f(θ)`。

★★**捩れの条件は落とせない**(落とすと偽になる)。
`Gal(K^{ab}/K) ≅ 𝒪_K^× × Ẑ` の分解自体が `π` に依存しており、
たとえば `θ := Art_{π'}(π')` は `ρ_{f'}(θ) = 1` だが `ρ_f(θ) = π'/π ≠ 1` である。
★捩れ元に限ると `Ẑ` 成分が捩れを持たないことから `θ` は慣性部分に入り、
慣性の上では Artin 写像は素元の選び方に依らない(これが Corollary 4.9 の内容)。
★★**したがってこの仮定は古典的に正しい。通すための偽の仮定ではない。**

★埋めるための在庫:`reciprocity_eq_of_intertwiner` /
`reciprocity_eq_of_thetaSet`(`LubinTateReciprocityIndependence.lean`、抽象核)、
`lubinTateClosure_sup_unramifiedClosure_eq_of_uniformizer`
(`LubinTateUniformizerIndependence.lean`)、`DworkMultiplicative.lean`。
★足りないのは `θ ∈ Θ^{K̂}_{π,π'}`(完備不分岐拡大の上の Lubin-Tate 同型)の構成である。 -/
def ReciprocityDatumIndependenceOnTorsion (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {π' : 𝒪[K.carrier]} (hπmax' : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π'})
    (hπne0' : π' ≠ 0)
    (f' : PowerSeries 𝒪[K.carrier]) (hf0' : PowerSeries.coeff 0 f' = 0)
    (hf1' : PowerSeries.coeff 1 f' = π')
    (hf' : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f' = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (θ : K.closure ≃ₐ[K.carrier] K.closure),
    (AlgEquiv.restrictNormalHom (abelianClosure K) θ) ^ (p ^ m) = 1 →
      ((reciprocityUnits K hq hπmax' hπne0' f' hf0' hf1' hf' θ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])
        = ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf θ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])

open scoped Classical in
/-- ★★★★★**穴 ⇒ 壁(1 つの `S` について)**。

★`E` は `abelianGalEquivUnitsZHat L hq hπmax hπne0 f ...`(★**構成的な方**)を取る。
`π` と `f` は `L` の上で任意に 1 つ選ぶ(`exists_lubinTateSeries`)。

★★**証明の流れ**(すべて本ファイルの前の節):
1. `x = τ|_{L^{ab}}` と書く(`restrictNormalHom` の全射性)。
2. `(E x).1 = ρ_f(τ)` / `(E (Ψ_g x)).1 = ρ_f(Ψ_g τ)`(§7、`E` は `ρ` 由来)。
3. `Ψ_g x` は `x` と同じく `p^n` 捩れ(`abelianGalConj` は群同型だから)。
4. 穴 `h` で `ρ_f(Ψ_g τ)` を `ρ_{f^{σ_g}}(Ψ_g τ)` に取り替える。
5. §6(`Γ_F` への代入)で `σ_g(ρ_f(τ))` に落ちる。

★★**`Ẑ` 成分には 1 度も触れていない** —— `ArtinUnitEquivariance` が要求するのが
単数成分だけだからである。 -/
theorem artinUnitEquivariance_at (h : ReciprocityDatumIndependenceOnTorsion p)
    (F : PAdicLocalField p) (n : ℕ) (S : Subgroup F.absGal) [S.Normal]
    (hopen : IsOpen (S : Set F.absGal)) (g : F.absGal) :
    ∃ E : Gal(↥(abelianClosure (fixedFieldLocalField F S hopen))/(fixedFieldLocalField F S hopen).carrier)
        ≃* ((𝒪[(fixedFieldLocalField F S hopen).carrier])ˣ × ZHat),
      ∀ x, x ^ (p ^ n) = 1 →
        ((unitsToField (fixedFieldLocalField F S hopen)
              (E (abelianGalConj F S hopen g x)).1 :
            ((fixedFieldLocalField F S hopen).carrier)ˣ)
          : (fixedFieldLocalField F S hopen).carrier)
          = fixedFieldAut S g
              ((unitsToField (fixedFieldLocalField F S hopen) (E x).1 :
                  ((fixedFieldLocalField F S hopen).carrier)ˣ)
                : (fixedFieldLocalField F S hopen).carrier) := by
  haveI := isAdicComplete_valuationRing (fixedFieldLocalField F S hopen)
  haveI := valuationRing_isDVR (fixedFieldLocalField F S hopen)
  obtain ⟨ϖ, hϖirr⟩ := IsDiscreteValuationRing.exists_irreducible
    𝒪[(fixedFieldLocalField F S hopen).carrier]
  have hπmax : IsLocalRing.maximalIdeal 𝒪[(fixedFieldLocalField F S hopen).carrier]
      = Ideal.span {ϖ} := (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖirr
  have hπne0 : ϖ ≠ 0 := hϖirr.ne_zero
  have hq : Fintype.card 𝓀[(fixedFieldLocalField F S hopen).carrier]
      = p ^ (absoluteInertiaDegree (fixedFieldLocalField F S hopen)) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow (fixedFieldLocalField F S hopen)
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries
    (A := 𝒪[(fixedFieldLocalField F S hopen).carrier]) hq hπmax
  refine ⟨abelianGalEquivUnitsZHat (fixedFieldLocalField F S hopen)
    hq hπmax hπne0 f hf0 hf1 hf, ?_⟩
  intro x hx
  obtain ⟨τ, rfl⟩ := AlgEquiv.restrictNormalHom_surjective
    (F := (fixedFieldLocalField F S hopen).carrier)
    (K₁ := (abelianClosure (fixedFieldLocalField F S hopen) : Type _))
    (fixedFieldLocalField F S hopen).closure x
  have htor : (AlgEquiv.restrictNormalHom (abelianClosure (fixedFieldLocalField F S hopen))
      (absGalConjCME F S hopen g τ)) ^ (p ^ n) = 1 := by
    rw [← abelianGalConj_restrictNormalHom, ← map_pow, hx, map_one]
  have hind := h (fixedFieldLocalField F S hopen) hq hπmax hπne0 f hf0 hf1 hf
    (maximalIdeal_eq_span_map (fixedFieldIntegerAut F S hopen g) hπmax)
    (fun hh => hπne0 ((fixedFieldIntegerAut F S hopen g).injective (by rw [hh, map_zero])))
    (PowerSeries.map ((fixedFieldIntegerAut F S hopen g :
        𝒪[(fixedFieldLocalField F S hopen).carrier] ≃+*
          𝒪[(fixedFieldLocalField F S hopen).carrier]) :
      𝒪[(fixedFieldLocalField F S hopen).carrier] →+*
        𝒪[(fixedFieldLocalField F S hopen).carrier]) f)
    (coeff_zero_map_twist _ hf0) (coeff_one_map_twist _ hf1)
    (map_residue_map_twist (fixedFieldIntegerAut F S hopen g) hf)
    n (absGalConjCME F S hopen g τ) htor
  rw [abelianGalConj_restrictNormalHom, abelianGalEquivUnitsZHat_restrictNormalHom_fst,
    abelianGalEquivUnitsZHat_restrictNormalHom_fst, unitsToField_coe, unitsToField_coe,
    ← hind, reciprocityUnits_absGalConjCME F S hopen hq g hπmax hπne0 f hf0 hf1 hf τ]
  exact coe_fixedFieldIntegerAut F S hopen g _

def artinUnitEquivariance_of_reciprocityDatumIndependenceOnTorsion.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- ★★★★★★**穴 ⇒ 壁**。

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

★`ArtinUnitEquivariance`(`ArtinEquivarianceProof.lean` の「残った壁」)が、
`ReciprocityDatumIndependenceOnTorsion` ただ 1 本から出る。

★★**無条件では出ていない。**残っているのは穴 1 本だけである。
★穴が埋まれば、`artinUnitEquivariance_iff_cyclotomeConj` と
`cyclotomicCharacter_recoverable_of_artinUnitEquivariance`(どちらも在庫)を
経由して、上の原文がそのまま出る。 -/
theorem artinUnitEquivariance_of_reciprocityDatumIndependenceOnTorsion
    (h : ReciprocityDatumIndependenceOnTorsion p) : ArtinUnitEquivariance p := by
  intro F n S hS hopen _hmu g
  exact @artinUnitEquivariance_at p _ h F n S hS hopen g

end ArtinReduction

end ABC3.Found.PGC
