import ABC3.Found.PGC.PowerSeriesAevalComm
import ABC3.Found.PGC.AdjoinIntegers

/-!
# 体の包含 `K.carrier⟮y⟯≤K.carrier⟮x⟯` に沿った冪級数評価の橋渡し(`sorry` 無し)

`LubinTateActionEquivariance.lean`はGalois同変性(`σ(a·x)=a·σ(x)`)を、
「`adjoinIntegers K x`と`adjoinIntegers K (σx)`という異なる2つの環を
橋渡しする(cross-point instance bridging)」を**回避する**ことで
確立した——`σ(x)`が`K.carrier⟮x⟯`自身に留まるという事実を使い、
評価を常に`x`自身の座標系だけで行った。

節目(5)(射影極限`Gal(L_π/K)≅𝒪_K^×`)へ向けては、この回避が使えない
場面が出る——`L_n≤L_{n+1}`という**体の包含**(自己同型ではない)に
沿って、`y`(`L_n`の生成元)自身の座標系での評価と、`x`(`L_{n+1}`の
生成元)の座標系で`y`を1点として評価したものが、同じ`K.closure`の
値を与えることを実際に示す必要がある——これが正面から向き合う必要の
ある cross-point bridging そのもの。

## 鍵となる一般化

`algHom_aeval_powerSeries_comm`(`PowerSeriesAevalComm.lean`)は
`σ:S→ₐ[A]S`という**自己**準同型にしか使えない設計だった
(`LubinTateActionEquivariance.lean`が自己同型限定の技法しか必要と
しなかったため)。証明そのものは`S`が2つに分かれても一字一句同じ
なので、`algHom_aeval_powerSeries_comm'`(始域`S₁`・終域`S₂`が
異なってよい版、同ファイルに追加済み)を`adjoinIntegers K y→ₐ[𝒪_K]
adjoinIntegers K x`という**埋め込み**(`IntermediateField.inclusion`
の制限)に適用するだけで、cross-point bridging が実際に組み上がる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **体の包含は`spectralNorm`を保つ**——`IntermediateField.inclusion`
は同じ`K.closure`の元をそのまま指す(`IntermediateField.coe_inclusion`)
だけなので、周囲のノルムをそのまま継承する。 -/
theorem norm_inclusion_eq {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (z : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    ‖IntermediateField.inclusion hle z‖ = ‖z‖ := by
  show spectralNorm K.carrier K.closure (↑(IntermediateField.inclusion hle z) : K.closure) =
    spectralNorm K.carrier K.closure (↑z : K.closure)
  rw [IntermediateField.coe_inclusion]

/-- **`adjoinIntegers K y`から`adjoinIntegers K x`への制限**——
ノルムを保つ(上の補題)ので、ノルム`≤1`の部分環をそのまま
`adjoinIntegers K x`(同じくノルム`≤1`の部分環)へ写す。 -/
noncomputable def adjoinIntegersInclusion {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    adjoinIntegers K y →+* adjoinIntegers K x where
  toFun z := ⟨IntermediateField.inclusion hle (z : IntermediateField.adjoin K.carrier ({y} : Set K.closure)),
    by
      show ‖IntermediateField.inclusion hle (z : IntermediateField.adjoin K.carrier ({y} : Set K.closure))‖ ≤ 1
      rw [norm_inclusion_eq]
      exact z.2⟩
  map_one' := Subtype.ext (map_one (IntermediateField.inclusion hle))
  map_mul' a b := Subtype.ext (map_mul (IntermediateField.inclusion hle) _ _)
  map_zero' := Subtype.ext (map_zero (IntermediateField.inclusion hle))
  map_add' a b := Subtype.ext (map_add (IntermediateField.inclusion hle) _ _)

/-- `adjoinIntegersInclusion`の座標——`K.closure`へ戻せば元の値と
一致する。 -/
theorem coe_adjoinIntegersInclusion {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (z : adjoinIntegers K y) :
    (↑(↑(adjoinIntegersInclusion K x y hle z) : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
        K.closure) =
      (↑(↑z : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) : K.closure) := by
  show (↑(IntermediateField.inclusion hle (z : IntermediateField.adjoin K.carrier ({y} : Set K.closure))) :
      K.closure) = _
  rw [IntermediateField.coe_inclusion]

/-- **`adjoinIntegersInclusion`は`𝒪_K`-代数準同型**——`IntermediateField.
inclusion`自身が`→ₐ[K.carrier]`であることの制限。`commutes'`は`K.closure`
へ戻せばどちらも`algebraMap K.carrier K.closure`そのものになることから
`rfl`で閉じる。 -/
noncomputable def adjoinIntegersInclusionAlgHom {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    adjoinIntegers K y →ₐ[𝒪[K.carrier]] adjoinIntegers K x where
  toFun := adjoinIntegersInclusion K x y hle
  map_one' := (adjoinIntegersInclusion K x y hle).map_one
  map_mul' := (adjoinIntegersInclusion K x y hle).map_mul
  map_zero' := (adjoinIntegersInclusion K x y hle).map_zero
  map_add' := (adjoinIntegersInclusion K x y hle).map_add
  commutes' c := by
    apply Subtype.ext; apply Subtype.ext
    rw [coe_adjoinIntegersInclusion]
    show algebraMap K.carrier K.closure (c : K.carrier) = algebraMap K.carrier K.closure (c : K.carrier)
    rfl

/-- **`adjoinIntegersInclusionAlgHom`は連続**——両側とも`K.carrier`
上有限次元なので、`K.carrier`-線形写像は自動的に連続
(`LinearMap.continuous_of_finiteDimensional`)。 -/
theorem continuous_adjoinIntegersInclusionAlgHom {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Continuous (adjoinIntegersInclusionAlgHom K x y hle) := by
  apply Continuous.subtype_mk
  exact (LinearMap.continuous_of_finiteDimensional (IntermediateField.inclusion hle).toLinearMap).comp
    continuous_subtype_val

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**cross-point bridging そのもの**: `y`自身の座標系(`adjoinIntegers
K y`)で冪級数を評価したものと、`x`の座標系で`y`(を`adjoinIntegers
K x`へ埋め込んだ点)を評価したものは、`K.closure`の値として一致する。
`algHom_aeval_powerSeries_comm'`(`PowerSeriesAevalComm.lean`、始域・
終域が異なってよい一般版)を`adjoinIntegersInclusionAlgHom`に適用
するだけ——`LubinTateActionEquivariance.lean`が自己同型限定の技法で
回避していた橋渡しを、埋め込みについて実際に構築した。 -/
theorem lubinTateEvalAtPoint_inclusion_comm {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure)
    (hle : IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (w : adjoinIntegers K y) (hw : PowerSeries.HasEval w)
    (hw' : PowerSeries.HasEval (adjoinIntegersInclusionAlgHom K x y hle w))
    (G : PowerSeries (𝒪[K.carrier])) :
    (↑(↑(lubinTateEvalAtPoint K x (adjoinIntegersInclusionAlgHom K x y hle w) hw' G) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) =
    (↑(↑(lubinTateEvalAtPoint K y w hw G) : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
        K.closure) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  haveI := completeSpace_adjoinIntegers K y
  haveI := isLinearTopology_adjoinIntegers K y
  haveI := continuousSMul_adjoinIntegers K y
  have hkey := algHom_aeval_powerSeries_comm' (adjoinIntegersInclusionAlgHom K x y hle)
    (continuous_adjoinIntegersInclusionAlgHom K x y hle) G hw hw'
  unfold lubinTateEvalAtPoint
  rw [← hkey]
  show (↑(↑((adjoinIntegersInclusion K x y hle) (PowerSeries.aeval hw G)) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = _
  rw [coe_adjoinIntegersInclusion]

end ABC3.Found.PGC
