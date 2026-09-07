import ABC3.Found.PGC.Prop22FilteredHypothesis
import ABC3.Found.PGC.PadicLogIntegers

/-!
# `reciprocityUnits` の **α-同変性** —— 何が無条件で、何が偽か

持ち場は「抽象 `α : Γ_K ≃ₜ* Γ_{K'}` に沿った Artin 写像(`reciprocityUnits`)の移送」。
在庫にあったのは `reciprocityUnits_semilinear_conj`(`σ : K ≃ₐ K` の**自己**同型版、
`Found/PGC/ReciprocityLimitEquivariance.lean`)と `reciprocityUnits_eq_of_transport`
(データ非依存性、`Found/PGC/ReciprocityDatumIndependence.lean`)だけで、
★**抽象 `α` 版は無かった**(測り方は下の「在庫調査」)。

## ★★★0. 結論を先に —— 3 行

1. ★**抽象核**: 全射 `φ : G ↠ H`・`φ' : G' ↠ H'` と同型 `α : G ≃* G'` について、
   「`α` に沿って `φ` を運ぶ `e : H ≃* H'` が在る」⟺「`Subgroup.map α (ker φ) = ker φ'`」
   (`exists_intertwines_iff`)。★**分岐・付値・Galois・Lubin-Tate が 1 語も出ない。**
2. ★★**したがって「`Art` の α-同変性」は仮説なしでは成り立たない。**
   ★`α = id` に取っても、両側の Lubin-Tate データが違えば偽になる ——
   `artinKerTransport_refl_iff_lubinTateClosure_eq` が
   ★★**「運べる ⟺ `K_π = K_{π′}`」を機械検査で示す**。そして `K_π = K_{π′}` は
   **偽**である(`Found/PGC/LubinTateUniformizerIndependence.lean` 冒頭の反例:
   `K = ℚ_p`(`p` 奇)、`π = p`、`π′ = -p`)。
3. ★★★★**核の対応(`ArtinKerTransport`)を仮説に置くと、そこから先は仮定ゼロで
   `𝒪_K ≃+ 𝒪_{K'}` まで一気に出る**(`integersTransport`)——
   これは pGC Proposition 2.2 の**第一段**そのものである。

## ★★1. 在庫調査(コマンドを残す)

```
grep -nE "ABC3\.Found\.PGC\.[a-zA-Z_]*(reciprocityUnits|artinMap|padicLog)[a-zA-Z_]*" \
  .cache/decl-index.txt
grep -nE "	(MulEquiv\.subgroupMap|MulEquiv\.subgroupCongr|QuotientGroup\.congr|
  QuotientGroup\.quotientKerEquivOfSurjective|AddEquiv\.toMultiplicative)	" \
  .cache/mathlib-index.txt
grep -nE "restrictNormalHom" .cache/mathlib-index.txt
```

見つかったもの(★「部品で引く」が効いた例):

| 要るもの | 在庫 |
|---|---|
| 第一同型定理 + 核の付け替え | `QuotientGroup.quotientKerEquivOfSurjective` / `QuotientGroup.congr`(mathlib) |
| 部分群への制限・部分群の付け替え | `MulEquiv.subgroupMap` / `MulEquiv.subgroupCongr`(mathlib) |
| `ker(restrictNormalHom E) = E.fixingSubgroup` | `IntermediateField.restrictNormalHom_ker`(mathlib) |
| `Art(Γ^n_K) = U^n_K`(仮定ゼロ) | `RamificationImageStage.lean` |
| `Art` の全射性 | `RamificationImageUnits.lean::reciprocityUnits_surjective` |
| ★`U^{r e_K}_K ≃* Multiplicative 𝒪_K`(**同型の形**) | `PadicLogIntegers.lean::padicLogPrincipalUnitsEquiv` |
| `Multiplicative` ↔ `Additive` の橋 | `AddEquiv.toMultiplicative`(mathlib) |

★★`padicLogPrincipalUnitsEquiv` が**集合の等式ではなく群同型**で在ったことが決定的だった。
持ち場の表が名指ししていた `smul_padicLog_image_ramificationFiltration_eq_integers`
(集合の等式)ではなく、★**その 1 つ手前の同型**を使っている。

## ★★2. 抽象核(§1)—— 型クラスは `Group` 4 つだけ

* `Intertwines φ φ' α e := ∀ g, e (φ g) = φ' (α g)`
* `ker_map_eq_of_intertwines` —— 絡めば核は対応する(★**全射性は要らない**)
* `transportOfKerEq` —— 核が対応すれば運べる
* `exists_intertwines_iff` —— ★両者は**同値**
* `intertwines_unique` —— 運び方は一意
* `map_map_of_intertwines` —— ★**部分群の族はそのまま運ばれる**(核の仮説すら不要)
* `intertwines_refl_conj` —— ★**標的が可換なら内部自己同型に沿った運びは恒等**

★分岐・付値・Galois・Lubin-Tate・`ℚ_p` は §1 に **1 語も出てこない**。

## ★★3. 具体層(§2〜§5)

`ArtinFilteredDatum K` は「`Γ_K ↠ 𝒪_K^×` で `Γ^n_K ↦ U^n_K` を満たす全射」+ 素元。
★**型に Lubin-Tate は出てこない**。`nonempty_artinFilteredDatum`(★**仮定ゼロ**)が
Lubin-Tate の選択(`π`・`f`)を `∃` の内側へ閉じ込める(`exists_artinMap` と同じ作法)。

仮説 `ArtinKerTransport D D' α := Subgroup.map α (ker Art_K) = ker Art_{K'}` の下で:

| 出るもの | 宣言 |
|---|---|
| `𝒪_K^× ≃* 𝒪_{K'}^×` と α-同変性 | `unitsTransport` / `unitsTransport_art` |
| ★`U^n_K ↦ U^n_{K'}`(**すべての `n`**) | `map_principalUnits_unitsTransport` |
| `(𝒪_K/π^n)^× ≅ (𝒪_{K'}/π′^n)^×` | `unitsQuotientTransport` |
| ★★★`𝒪_K ≃+ 𝒪_{K'}` | `integersTransport` |

★`U^n` の対応には **α が濾過つき**であることを使う(`FilteredGroup.Iso.map_Gv`)。
★`𝒪_K ≃+ 𝒪_{K'}` には `e_K = e_{K'}` は**要らない** —— 段の番号に公倍数
`2·e_K·e_{K'}` を使う(`integersTransport` の docstring)。

## ★★★4. 退化の自己検査(隠さずに書く)

1. ★**非空虚**: 恒等(`artinKerTransport_refl`)と**内部自己同型**
   (`artinKerTransport_inner`)では仮説が**無条件に**成り立つ。
   ★内部自己同型では `𝒪_K^×` が可換なので、運びは**恒等写像**である
   (`unitsTransport_inner`)——「非空虚だが内容が薄い場合」がここに見えている。
2. ★★**自明ではない**: `artinKerTransport_refl_iff_lubinTateClosure_eq` が示すとおり、
   ★`α = id` でもデータが違えば偽になる。★**仮説は本物である。**
3. ★★**`𝒪_K ≃+ 𝒪_{K'}` という結論だけを見ると弱い** —— 加法群としては
   `𝒪_K ≅ ℤ_p^{[K:ℚ_p]}` なので、次数が等しければ **α と無関係に**同型である。
   ★本ファイルの内容は「同型が在る」ことではなく
   ★**「`α` から**標準的に**作れて、しかも `U^n` の対応を伴う」**ことである
   (`map_principalUnits_unitsTransport` が全 `n` で成り立つのが実質的な中身)。
   ★`n = 1` を入れると `k_K^× ≅ k_{K'}^×`(`residueUnitsTransport`)——
   これは `α` に本当の制約を課す(剰余体の位数が一致する)。
4. ★★**まだ Proposition 2.2 ではない**: `IntKbarTransportFiltered` が要求するのは
   ★`𝒪_{K̄} ≃+ 𝒪_{K̄'}` の **Γ-同変**同型であり、本ファイルが作るのは
   **底体の** `𝒪_K ≃+ 𝒪_{K'}`(Γ は自明に作用)である。★**渡せる形にはなっていない。**
   次の節点は「開部分群 ↔ 有限次部分拡大」を通して `𝒪_L ≃+ 𝒪_{L'}` を整合的に作り、
   `𝒪_{K̄} = colim_L 𝒪_L` を取ること(下の「次の節点」)。

## ★★5. 原典より短い道(名指し)

★原典 Proposition 2.2 の証明は「有限次拡大 `L/K` に降りて Prop 2.1 を使い、
上付き→下付き→上付きと番号付けを往復する」。★本ファイルは
★**番号付けの往復を 1 度も行わない** —— `Art(Γ^n) = U^n` を上付きのまま使い、
`log` で `𝒪` に落とすだけである。
★さらに `e_K = e_{K'}` を示す代わりに ★**段の番号に公倍数を使う**ことで、
「絶対分岐指数が `α` から復元できるか」という(未解決の)問いを**回避した**。

## ★6. 次の節点(名指し)

1. ★`α` が開部分群 `H ≤ Γ_K` を開部分群 `α(H) ≤ Γ_{K'}` に写すことから
   有限次部分拡大 `L ↔ L'` の対応を作り、各段で `ArtinFilteredDatum L` / `ArtinFilteredDatum L'` と
   `ArtinKerTransport` を**整合的に**取る(★整合性が非自明な部分)。
2. ★`𝒪_{K̄} = colim_L 𝒪_L` の余極限を取って Γ-同変性を出す。
3. ★`ArtinKerTransport` を**落とせるか**: `ker Art_K ∩ I_K` は
   `Gal(K̄/K^{ab}) ∩ I_K`(=閉包した交換子群 ∩ 惰性群)に等しいはずで、これは
   **群論的に標準**である。★これが言えれば「惰性群に制限した Artin 写像」の
   α-同変性は**仮説なしで**出る。要るのは局所 Kronecker-Weber
   (`LocalClassFieldTheory.lean`、木にある)と `K_π ⊔ K^ur = K^{ab}`。
   ★**本シフトでは測っていない**(着手していない)。

## ★衝突検査(#158)

★実名前空間に積んだところ `lake build ABC3.Found` が

```
error: ABC3/Found.lean:1:0: import ABC3.Found.PGC.ReciprocityAlphaTransport failed,
environment already contains 'ABC3.Found.PGC.ArtinDatum.mk.noConfusion'
from ABC3.Found.PGC.Section3RealParameters
```

を返した。★`Found/PGC/Section3RealParameters.lean:231` の `ArtinDatum`
(`K^× →* Gal(K^ab/K)` の側、分岐濾過の条件を持たない)と**別物**なので、
本ファイルの構造体は `ArtinFilteredDatum` に改名した。
★**ファイル単体の `leanfile.mjs` では出ない**(相手を import していないから)。

## 逸脱の記録

1. `Skeleton/**` は 1 行も触っていない。`Found/PGC/LubinTate*.lean` も読むだけ(D28)。
2. 原典は `Art` を `Γ^{ab}_K` の言葉で述べるが、木の `reciprocityUnits` は
   `𝒪_K^×` 成分だけを取る(`RamificationImageStage.lean` の「逸脱の記録 2」を踏襲)。
3. 原典 Proposition 2.2 は `v = r·e_K`(`r ≥ 2`)を使うが、本ファイルは
   `v = 2·e_K·e_{K'}` を使う(上記 §5)。★原典の条件 `r ≥ 2` は両側で満たしている。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open IsLocalRing
open scoped NNReal Valued

/-! ## §1 ★★★抽象核 —— **純群論**。分岐・付値・Galois・Lubin-Tate は 1 語も出ない -/

section TransportCore

variable {G G' H H' : Type*} [Group G] [Group G'] [Group H] [Group H']

/-- `e` が `φ` と `φ'` を `α` に沿って絡める、という関係。 -/
def Intertwines (φ : G →* H) (φ' : G' →* H') (α : G ≃* G') (e : H ≃* H') : Prop :=
  ∀ g : G, e (φ g) = φ' (α g)

/-- ★**絡む `e` があれば、`α` は核を核へちょうど写す**(全射性は要らない)。 -/
theorem ker_map_eq_of_intertwines {φ : G →* H} {φ' : G' →* H'} {α : G ≃* G'} {e : H ≃* H'}
    (he : Intertwines φ φ' α e) : Subgroup.map (α : G →* G') φ.ker = φ'.ker := by
  ext y
  simp only [Subgroup.mem_map, MonoidHom.mem_ker, MonoidHom.coe_coe]
  constructor
  · rintro ⟨g, hg, rfl⟩
    rw [← he g, hg, map_one]
  · intro hy
    refine ⟨α.symm y, ?_, by simp⟩
    have hkey := he (α.symm y)
    rw [MulEquiv.apply_symm_apply, hy] at hkey
    exact e.injective (by rw [hkey, map_one])

/-- ★★**核が対応するなら、全射準同型は `α` に沿って運べる**。
第一同型定理 2 回 + `QuotientGroup.congr`。 -/
noncomputable def transportOfKerEq (φ : G →* H) (hφ : Function.Surjective φ)
    (φ' : G' →* H') (hφ' : Function.Surjective φ') (α : G ≃* G')
    (hker : Subgroup.map (α : G →* G') φ.ker = φ'.ker) : H ≃* H' :=
  ((QuotientGroup.quotientKerEquivOfSurjective φ hφ).symm.trans
      (QuotientGroup.congr φ.ker φ'.ker α hker)).trans
    (QuotientGroup.quotientKerEquivOfSurjective φ' hφ')

theorem transportOfKerEq_intertwines (φ : G →* H) (hφ : Function.Surjective φ)
    (φ' : G' →* H') (hφ' : Function.Surjective φ') (α : G ≃* G')
    (hker : Subgroup.map (α : G →* G') φ.ker = φ'.ker) :
    Intertwines φ φ' α (transportOfKerEq φ hφ φ' hφ' α hker) := by
  intro g
  have h1 : (QuotientGroup.quotientKerEquivOfSurjective φ hφ).symm (φ g)
      = (QuotientGroup.mk g : G ⧸ φ.ker) := by
    rw [MulEquiv.symm_apply_eq]
    exact (QuotientGroup.kerLift_mk φ g).symm
  show (QuotientGroup.quotientKerEquivOfSurjective φ' hφ')
      (QuotientGroup.congr φ.ker φ'.ker α hker
        ((QuotientGroup.quotientKerEquivOfSurjective φ hφ).symm (φ g))) = _
  rw [h1, QuotientGroup.congr_mk]
  exact QuotientGroup.kerLift_mk φ' (α g)

/-- ★絡む `e` は(`φ` が全射なら)**一意**。 -/
theorem intertwines_unique {φ : G →* H} (hφ : Function.Surjective φ) {φ' : G' →* H'}
    {α : G ≃* G'} {e₁ e₂ : H ≃* H'} (h₁ : Intertwines φ φ' α e₁)
    (h₂ : Intertwines φ φ' α e₂) : e₁ = e₂ := by
  ext h
  obtain ⟨g, rfl⟩ := hφ h
  rw [h₁ g, h₂ g]

/-- ★★★**抽象核の到達点** —— 「`α` に沿って運べる」⟺「`α` が核を対応させる」。 -/
theorem exists_intertwines_iff (φ : G →* H) (hφ : Function.Surjective φ)
    (φ' : G' →* H') (hφ' : Function.Surjective φ') (α : G ≃* G') :
    (∃ e : H ≃* H', Intertwines φ φ' α e) ↔ Subgroup.map (α : G →* G') φ.ker = φ'.ker :=
  ⟨fun ⟨_, he⟩ => ker_map_eq_of_intertwines he,
    fun h => ⟨transportOfKerEq φ hφ φ' hφ' α h, transportOfKerEq_intertwines φ hφ φ' hφ' α h⟩⟩

/-- ★★**部分群の族はそのまま運ばれる**(核の仮説すら要らない —— 絡めば出る)。 -/
theorem map_map_of_intertwines {φ : G →* H} {φ' : G' →* H'} {α : G ≃* G'} {e : H ≃* H'}
    (he : Intertwines φ φ' α e) (S : Subgroup G) :
    Subgroup.map (e : H →* H') (Subgroup.map φ S)
      = Subgroup.map φ' (Subgroup.map (α : G →* G') S) := by
  rw [Subgroup.map_map, Subgroup.map_map]
  congr 1
  exact MonoidHom.ext he

/-- ★★**標的が可換なら、内部自己同型に沿った運びは恒等**。
`φ (c g c⁻¹) = φ c · φ g · (φ c)⁻¹ = φ g`。 -/
theorem intertwines_refl_conj {H : Type*} [CommGroup H] (φ : G →* H) (α : G ≃* G) (c : G)
    (hα : ∀ g, α g = c * g * c⁻¹) : Intertwines φ φ α (MulEquiv.refl H) := by
  intro g
  simp only [MulEquiv.refl_apply, hα, map_mul, map_inv]
  rw [mul_comm (φ c) (φ g), mul_assoc, mul_inv_cancel, mul_one]

/-- ★恒等は常に絡む。 -/
theorem intertwines_refl (φ : G →* H) : Intertwines φ φ (MulEquiv.refl G) (MulEquiv.refl H) :=
  fun _ => rfl

end TransportCore

/-! ## §2 具体層 —— **Artin データ**(Lubin-Tate の選択を 1 つの構造に閉じ込める) -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★**`Γ_K` の上の Artin データ**。

`Γ_K` から `𝒪_K^×` への**全射**準同型で、**上付き分岐濾過を主単数濾過へ写す**もの。
★型に Lubin-Tate は 1 語も出てこない(`unif` は素元というだけ)。 -/
structure ArtinFilteredDatum (K : PAdicLocalField p) where
  /-- 相互写像 `Art_K : Γ_K ↠ 𝒪_K^×`。 -/
  art : K.absGal →* (𝒪[K.carrier])ˣ
  /-- 全射(`RamificationImageUnits.lean::reciprocityUnits_surjective`)。 -/
  art_surjective : Function.Surjective art
  /-- 素元。 -/
  unif : 𝒪[K.carrier]
  /-- `unif` は本当に素元。 -/
  unif_max : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {unif}
  /-- 素元は 0 でない。 -/
  unif_ne_zero : unif ≠ 0
  /-- ★`Art(Γ_K^n) = U^n_K`(`RamificationImageStage.lean`、仮定ゼロ)。 -/
  map_Gv : ∀ n : ℕ, Subgroup.map art ((ramificationFiltration p).Gv K (n : ℝ))
    = principalUnits K unif n

def ArtinFilteredDatum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★★**Artin データは無条件に存在する**(`exists_artinMap` と同じ作法で
Lubin-Tate の選択を `∃` の内側へ閉じ込める)。 -/
theorem nonempty_artinFilteredDatum (K : PAdicLocalField p) : Nonempty (ArtinFilteredDatum K) := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖirr
  have hπne0 : ϖ ≠ 0 := hϖirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  exact ⟨{ art := reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf
           art_surjective := reciprocityUnits_surjective K hq hπmax hπne0 f hf0 hf1 hf
           unif := ϖ
           unif_max := hπmax
           unif_ne_zero := hπne0
           map_Gv := fun n =>
             map_ramificationFiltration_reciprocityUnits_eq_principalUnits
               K hq hπmax hπne0 f hf0 hf1 hf n }⟩

def nonempty_artinFilteredDatum.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★**退化していない**: `Art(I_K) = 𝒪_K^×`(`n = 0` を入れるだけ)。 -/
theorem ArtinFilteredDatum.map_absInertia {K : PAdicLocalField p} (D : ArtinFilteredDatum K) :
    Subgroup.map D.art (absInertia K) = ⊤ := by
  have h := D.map_Gv 0
  rwa [Nat.cast_zero, ramificationFiltration_Gv_zero K, principalUnits_zero_eq_top K] at h

/-! ## §3 α に沿った移送 -/

/-- ★★**仮説** —— `α` が 2 つの Artin データの**核**を対応させる。 -/
def ArtinKerTransport {K K' : PAdicLocalField p} (D : ArtinFilteredDatum K) (D' : ArtinFilteredDatum K')
    (α : ContinuousMulEquiv K.absGal K'.absGal) : Prop :=
  Subgroup.map (α.toMulEquiv : K.absGal →* K'.absGal) D.art.ker = D'.art.ker

def ArtinKerTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**単数群の移送** `𝒪_K^× ≃* 𝒪_{K'}^×`(抽象核の代入)。 -/
noncomputable def unitsTransport {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K}
    {D' : ArtinFilteredDatum K'} {α : ContinuousMulEquiv K.absGal K'.absGal}
    (h : ArtinKerTransport D D' α) : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ :=
  transportOfKerEq D.art D.art_surjective D'.art D'.art_surjective α.toMulEquiv h

theorem unitsTransport_intertwines {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K}
    {D' : ArtinFilteredDatum K'} {α : ContinuousMulEquiv K.absGal K'.absGal}
    (h : ArtinKerTransport D D' α) :
    Intertwines D.art D'.art α.toMulEquiv (unitsTransport h) :=
  transportOfKerEq_intertwines _ _ _ _ _ h

/-- ★★**α-同変性そのもの**: `e (Art_K g) = Art_{K'} (α g)`。 -/
theorem unitsTransport_art {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K} {D' : ArtinFilteredDatum K'}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : ArtinKerTransport D D' α) (g : K.absGal) :
    unitsTransport h (D.art g) = D'.art (α g) :=
  unitsTransport_intertwines h g

/-- ★★★**核の対応は、α-同変な単数群同型の存在と同値**。

★これが本ファイルの中心である。「`Art` を `α` に沿って運べるか」は
**核が対応するかどうか、それだけ**で決まる。 -/
theorem exists_unitsEquiv_iff_artinKerTransport {K K' : PAdicLocalField p} (D : ArtinFilteredDatum K)
    (D' : ArtinFilteredDatum K') (α : ContinuousMulEquiv K.absGal K'.absGal) :
    (∃ e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ,
        ∀ g : K.absGal, e (D.art g) = D'.art (α g))
      ↔ ArtinKerTransport D D' α :=
  exists_intertwines_iff D.art D.art_surjective D'.art D'.art_surjective α.toMulEquiv

def exists_unitsEquiv_iff_artinKerTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★運び方は一意。 -/
theorem unitsTransport_unique {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K} {D' : ArtinFilteredDatum K'}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : ArtinKerTransport D D' α)
    (e : (𝒪[K.carrier])ˣ ≃* (𝒪[K'.carrier])ˣ)
    (he : ∀ g : K.absGal, e (D.art g) = D'.art (α g)) : e = unitsTransport h :=
  intertwines_unique D.art_surjective he (unitsTransport_intertwines h)

/-! ### ★★濾過つき `α` なら、主単数濾過もそのまま運ばれる -/

/-- ★★★**`U^n_K` は `U^n_{K'}` へ写る**。★核の仮説と、`α` が濾過つきであることの
2 つだけを使う。 -/
theorem map_principalUnits_unitsTransport {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K}
    {D' : ArtinFilteredDatum K'} {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}
    (h : ArtinKerTransport D D' (filteredIsoEquiv α)) (n : ℕ) :
    Subgroup.map (unitsTransport h : (𝒪[K.carrier])ˣ →* (𝒪[K'.carrier])ˣ)
        (principalUnits K D.unif n) = principalUnits K' D'.unif n := by
  have hmap := map_map_of_intertwines (unitsTransport_intertwines h)
    ((ramificationFiltration p).Gv K ((n : ℕ) : ℝ))
  rw [D.map_Gv n] at hmap
  rw [hmap, show Subgroup.map ((filteredIsoEquiv α).toMulEquiv : K.absGal →* K'.absGal)
      ((ramificationFiltration p).Gv K ((n : ℕ) : ℝ))
      = (ramificationFiltration p).Gv K' ((n : ℕ) : ℝ) from α.map_Gv ((n : ℕ) : ℝ)]
  exact D'.map_Gv n

def map_principalUnits_unitsTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**系**: `(𝒪_K/π^n)^× ≅ (𝒪_{K'}/π'^n)^×`。
主単数濾過が対応するので、剰余環の単数群がそのまま対応する。 -/
noncomputable def unitsQuotientTransport {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K}
    {D' : ArtinFilteredDatum K'} {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}
    (h : ArtinKerTransport D D' (filteredIsoEquiv α)) (n : ℕ) (hn : 1 ≤ n) :
    (𝒪[K.carrier] ⧸ Ideal.span ({D.unif ^ n} : Set 𝒪[K.carrier]))ˣ
      ≃* (𝒪[K'.carrier] ⧸ Ideal.span ({D'.unif ^ n} : Set 𝒪[K'.carrier]))ˣ :=
  ((principalUnitsQuotientEquiv K D.unif_max n hn).symm.trans
      (QuotientGroup.congr (principalUnits K D.unif n) (principalUnits K' D'.unif n)
        (unitsTransport h) (map_principalUnits_unitsTransport h n))).trans
    (principalUnitsQuotientEquiv K' D'.unif_max n hn)

/-- ★★**退化していないことの witness**(`n = 1`)—— 剰余体の単数群が対応する:
`k_K^× ≅ k_{K'}^×`。★これは `α` に本当の制約を課す(剰余体の位数が一致する)。 -/
noncomputable def residueUnitsTransport {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K}
    {D' : ArtinFilteredDatum K'} {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}
    (h : ArtinKerTransport D D' (filteredIsoEquiv α)) :
    (𝒪[K.carrier] ⧸ Ideal.span ({D.unif ^ 1} : Set 𝒪[K.carrier]))ˣ
      ≃* (𝒪[K'.carrier] ⧸ Ideal.span ({D'.unif ^ 1} : Set 𝒪[K'.carrier]))ˣ :=
  unitsQuotientTransport h 1 le_rfl

/-! ## §4 ★★非空虚性 —— 仮説が**無条件に**成り立つ `α` -/

/-- ★恒等では無条件。 -/
theorem artinKerTransport_refl {K : PAdicLocalField p} (D : ArtinFilteredDatum K) :
    ArtinKerTransport D D (ContinuousMulEquiv.refl K.absGal) :=
  ker_map_eq_of_intertwines (intertwines_refl D.art)

/-- ★★**内部自己同型では無条件**。標的 `𝒪_K^×` が**可換**だから
`Art(c g c⁻¹) = Art g` になる —— 核が正規であること以上の内容がある。 -/
theorem artinKerTransport_inner {K : PAdicLocalField p} (D : ArtinFilteredDatum K) (c : K.absGal) :
    ArtinKerTransport D D (innerAbsGalEquiv K c) :=
  ker_map_eq_of_intertwines
    (intertwines_refl_conj D.art (innerAbsGalEquiv K c).toMulEquiv c (fun _ => rfl))

/-- ★★内部自己同型に沿った移送は**恒等写像**である。 -/
theorem unitsTransport_inner {K : PAdicLocalField p} (D : ArtinFilteredDatum K) (c : K.absGal) :
    unitsTransport (artinKerTransport_inner D c) = MulEquiv.refl (𝒪[K.carrier])ˣ :=
  (unitsTransport_unique _ _
    (intertwines_refl_conj D.art (innerAbsGalEquiv K c).toMulEquiv c (fun _ => rfl))).symm

/-- ★濾過つきの形で: 恒等。 -/
theorem artinKerTransportFiltered_refl {K : PAdicLocalField p} (D : ArtinFilteredDatum K) :
    ArtinKerTransport D D (filteredIsoEquiv (filteredIsoRefl (pgcFilteredGroup K))) :=
  artinKerTransport_refl D

/-- ★濾過つきの形で: 内部自己同型。 -/
theorem artinKerTransportFiltered_inner {K : PAdicLocalField p} (D : ArtinFilteredDatum K)
    (c : K.absGal) :
    ArtinKerTransport D D (filteredIsoEquiv (filteredIsoInner (pgcFilteredGroup K) c)) :=
  artinKerTransport_inner D c

/-! ## §5 ★★★★★`𝒪_K ≃+ 𝒪_{K'}` —— pGC Proposition 2.2 の第一段が α に沿って運べる -/

section IntegersTransport

/-- `e_L ≥ 1` なので `2 ≤ 2·e_L`。★2 つの体の `e` が**等しい必要はない**:
公倍数 `2·e_K·e_{K'}` を段の番号に使えばよい。 -/
theorem two_le_two_mul_absoluteRamificationIndex (L : PAdicLocalField p) :
    2 ≤ 2 * absoluteRamificationIndex L := by
  have h := absoluteRamificationIndex_ne_zero L
  omega

variable {K K' : PAdicLocalField p} {D : ArtinFilteredDatum K} {D' : ArtinFilteredDatum K'}
  {α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')}

/-- ★★★★★**`𝒪_K` と `𝒪_{K'}` は加法群として対応する**。

経路(すべて仮定ゼロの在庫の合成):

```
𝒪_K  ≃(log)  U^{2 e_K e_{K'}}_K  ≃(Art^{-1})  Γ_K^{2 e_K e_{K'}}
     ≃(α)    Γ_{K'}^{2 e_K e_{K'}}  ≃(Art)  U^{2 e_K e_{K'}}_{K'}  ≃(log)  𝒪_{K'}
```

★★段の番号に **`2·e_K·e_{K'}`(公倍数)** を使うので `e_K = e_{K'}` は要らない。
`padicLogPrincipalUnitsEquiv` が要求する `v = r·e` の形は、`K` 側では
`r = 2 e_{K'}`、`K'` 側では `r = 2 e_K` として満たされる(どちらも `≥ 2`)。

★使った仮説は **`ArtinKerTransport`(核の対応)1 つだけ**である。 -/
noncomputable def integersTransport (h : ArtinKerTransport D D' (filteredIsoEquiv α)) :
    𝒪[K.carrier] ≃+ 𝒪[K'.carrier] :=
  AddEquiv.toMultiplicative.symm
    ((padicLogPrincipalUnitsEquiv K D.unif_max D.unif_ne_zero
        (two_le_two_mul_absoluteRamificationIndex K')).symm.trans
      (((unitsTransport h).subgroupMap
          (principalUnits K D.unif
            (2 * absoluteRamificationIndex K' * absoluteRamificationIndex K))).trans
        ((MulEquiv.subgroupCongr (map_principalUnits_unitsTransport h _)).trans
          ((MulEquiv.subgroupCongr
              (congrArg (principalUnits K' D'.unif) (by ring))).trans
            (padicLogPrincipalUnitsEquiv K' D'.unif_max D'.unif_ne_zero
              (two_le_two_mul_absoluteRamificationIndex K))))))

def integersTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★消費しやすい形。 -/
theorem nonempty_integers_addEquiv_of_artinKerTransport
    (h : ArtinKerTransport D D' (filteredIsoEquiv α)) :
    Nonempty (𝒪[K.carrier] ≃+ 𝒪[K'.carrier]) := ⟨integersTransport h⟩

end IntegersTransport

/-! ## §6 ★★★**仮説は落とせない** —— `α = id` でも中身は `K_π = K_{π′}` である -/

theorem subgroup_map_mulEquivRefl {G : Type*} [Group G] (S : Subgroup G) :
    Subgroup.map ((MulEquiv.refl G : G ≃* G) : G →* G) S = S := by
  ext g
  simp

/-- ★★**`α = id` での仮説は「2 つの Artin データの核が一致する」ことに等しい**。 -/
theorem artinKerTransport_refl_iff {K : PAdicLocalField p} (D D' : ArtinFilteredDatum K) :
    ArtinKerTransport D D' (ContinuousMulEquiv.refl K.absGal) ↔ D.art.ker = D'.art.ker := by
  rw [ArtinKerTransport,
    show ((ContinuousMulEquiv.refl K.absGal).toMulEquiv : K.absGal →* K.absGal)
      = ((MulEquiv.refl K.absGal : K.absGal ≃* K.absGal) : K.absGal →* K.absGal) from rfl,
    subgroup_map_mulEquivRefl]

section LubinTateDatum

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

/-- Lubin-Tate のデータから作った Artin データ。 -/
noncomputable def lubinTateArtinFilteredDatum : ArtinFilteredDatum K where
  art := reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf
  art_surjective := reciprocityUnits_surjective K hq hπmax hπne0 f hf0 hf1 hf
  unif := π
  unif_max := hπmax
  unif_ne_zero := hπne0
  map_Gv n :=
    map_ramificationFiltration_reciprocityUnits_eq_principalUnits K hq hπmax hπne0 f hf0 hf1 hf n

@[simp] theorem lubinTateArtinFilteredDatum_art :
    (lubinTateArtinFilteredDatum K hq hπmax hπne0 f hf0 hf1 hf).art
      = reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf := rfl

variable [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]

/-- ★★**`Art_π` の核は `K_π` の固定部分群** `Gal(K̄/K_π)`。

`reciprocityUnits = Φ ∘ (K_π への制限)` で `Φ` は同型だから、核は制限の核に等しく、
それは mathlib の `IntermediateField.restrictNormalHom_ker` で固定部分群になる。 -/
theorem ker_reciprocityUnits :
    (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker
      = (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup := by
  rw [← IntermediateField.restrictNormalHom_ker
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
  ext σ
  simp only [MonoidHom.mem_ker]
  rw [← lubinTateClosureGalEquivUnits_restrictNormalHom K hq hπmax hπne0 f hf0 hf1 hf σ]
  exact ⟨fun h => (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf).injective
      (by rw [h, map_one]), fun h => by rw [h, map_one]⟩

def ker_reciprocityUnits.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

/-- ★★★★**本ファイルの一番大事な帰結** ——

同じ `K` の上の 2 つの Lubin-Tate データについて、
★**`α = id` で `Art` を運べる ⟺ `K_π = K_{π′}`**。

★★`K_π = K_{π′}` は **偽**である
(`Found/PGC/LubinTateUniformizerIndependence.lean` 冒頭の反例:
`K = ℚ_p`(`p` 奇素数)、`π = p`、`π′ = -p`。`K_p = ℚ_p(µ_{p^∞})` だが
`Art(-1)` は `µ_{p^∞}` 上で反転として働くので `K_{-p} ≠ K_p`)。

⇒ ★★★**「`α` に沿って `Art` が運べる」を仮説なしで主張することはできない。**
★恒等 `α`(濾過つき同型としては最も従順なもの)ですら、データの選び方次第で偽になる。
★★したがって `ArtinKerTransport` は**空虚でもなければ自明でもない**本物の仮説である。 -/
theorem artinKerTransport_refl_iff_lubinTateClosure_eq
    {π' : 𝒪[K.carrier]} (hπ'max : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π'})
    (hπ'ne0 : π' ≠ 0)
    (f' : PowerSeries 𝒪[K.carrier]) (hf'0 : PowerSeries.coeff 0 f' = 0)
    (hf'1 : PowerSeries.coeff 1 f' = π')
    (hf' : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f' = PowerSeries.X ^ (pp ^ ff))
    [Normal K.carrier (lubinTateClosure K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf')] :
    ArtinKerTransport (lubinTateArtinFilteredDatum K hq hπmax hπne0 f hf0 hf1 hf)
        (lubinTateArtinFilteredDatum K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf')
        (ContinuousMulEquiv.refl K.absGal)
      ↔ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
          = lubinTateClosure K hq hπ'max hπ'ne0 f' hf'0 hf'1 hf' := by
  haveI := isGalois_closure K
  rw [artinKerTransport_refl_iff, lubinTateArtinFilteredDatum_art, lubinTateArtinFilteredDatum_art,
    ker_reciprocityUnits, ker_reciprocityUnits]
  refine ⟨fun h => ?_, fun h => by rw [h]⟩
  rw [← InfiniteGalois.fixedField_fixingSubgroup
      (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf), h,
    InfiniteGalois.fixedField_fixingSubgroup]

def artinKerTransport_refl_iff_lubinTateClosure_eq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

end LubinTateDatum

#print axioms Intertwines
#print axioms transportOfKerEq
#print axioms exists_intertwines_iff
#print axioms map_map_of_intertwines
#print axioms intertwines_refl_conj
#print axioms nonempty_artinFilteredDatum
#print axioms exists_unitsEquiv_iff_artinKerTransport
#print axioms map_principalUnits_unitsTransport
#print axioms unitsQuotientTransport
#print axioms artinKerTransport_inner
#print axioms unitsTransport_inner
#print axioms artinKerTransport_refl
#print axioms residueUnitsTransport
#print axioms integersTransport
#print axioms nonempty_integers_addEquiv_of_artinKerTransport
#print axioms two_le_two_mul_absoluteRamificationIndex
#print axioms artinKerTransport_refl_iff
#print axioms ker_reciprocityUnits
#print axioms artinKerTransport_refl_iff_lubinTateClosure_eq
#print axioms lubinTateArtinFilteredDatum
#print axioms ArtinFilteredDatum.map_absInertia
#print axioms unitsTransport_art
#print axioms intertwines_unique
#print axioms ker_map_eq_of_intertwines
#print axioms transportOfKerEq_intertwines

end ABC3.Found.PGC
