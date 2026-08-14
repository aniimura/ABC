import ABC3.Check.PGC.DischargeFiresNothing

/-!
# [pGC] §1 の3定理は、現在の `Interface` の下で反証できるか — 横断監査

`Check/IUTchIII/Cor312Degenerate.lean` で、Corollary 3.12 が
**現在の `Interface` の下で偽**であることが分かった(反例を構成済み、`sorry` 無し)。
新しい失敗形を見つけたら既存の資産を横断的に再監査する、という規律に従い、
pGC §1 の3定理について同じことを調べる。

## 結果 — **3件とも反証できなかった**(0/3)

| 定理 | `Interface` | 反証 | 何が塞いでいるか |
|---|---|---|---|
| `cyclotomicCharacter_recoverable` | 取らない | **できない** | 共役の経路が**証明として閉じている**(下記1) |
| `residueCard_and_degree_recoverable` | `ResidueCardinality` | できない | K′=K の経路は自明に閉じる。残るのは α の構成のみ(下記2) |
| `inertia_recoverable` | `ResidueCardinality`・`SubgroupCorrespondence` | できない | 構成できる唯一の `SC` では**定理が証明できてしまう**(下記3) |

**これは「pGC の3定理が真である」という意味ではない。** 我々が反例を**構成できない**というだけである
(`LeanStatus.absent (searched)` と同じ規律——探した範囲を下に書く)。

## ★Corollary 3.12 との違い — 何が反証可能性を決めているか

表面的には「`Interface` が抽象的か具体的か」に見える。Cor 3.12 の `Interface` は
`Amb : Type` と `logVol : Set Amb → WithTop ℝ` という**自由形式の抽象データ**なので、
`Amb := Unit`・`logVol := (空集合なら +∞)` とその場で反例が書ける。pGC の `Interface` は
`PAdicLocalField p` や `Subgroup K.absGal` という**具体的な mathlib 対象**の上にある。

だが、より正確な言い方がある——**statement の量化の形**である:

```
Cor 3.12 :  ∀ D : PilotObjectData, (D だけで決まる結論)
pGC       :  ∀ RD SC, ∀ {K K'} (α : ContinuousMulEquiv K.absGal K'.absGal), (…)
```

Cor 3.12 の結論は **`Interface` のデータだけで決まる**。データが自由なら結論も自由に動かせる
——反証は `Interface` を書き換えるだけで済む。

pGC の結論は `Interface` の外の対象 **K・K′・α** にも量化している。反証するには
それらを**構成しなければならない**。しかも `refutation_reduces_to_alpha` が示すとおり、
反証に要る α は**その定理自身が存在しないと主張しているもの**である。

> **反証が効くのは、statement の falsity が `Interface` のデータだけで目撃できるときに限る。**
> 効かないことは「その定理が健全である」ことを意味しない——
> 我々が反例を**構成できない**というだけである。

この非対称性が G7(反証可能性の検査)を制度化するかどうかの判断材料になる。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## 我々が構成できる唯一の非自明な α — 共役 -/

/-- **共役による連続自己同型** `x ↦ h x h⁻¹`。

`RecoverableFromAbsGal` は「任意の連続同型 α について」と量化しているので、
反証にはまず α が要る。恒等射(`ContinuousMulEquiv.refl`)以外で我々が構成できるのは
これだけである——Krull 位相が位相群であること(`IsTopologicalGroup Gal(L/K)`、
mathlib `FieldTheory/KrullTopology.lean:143`)から連続性が出る。 -/
noncomputable def conjEquiv (K : PAdicLocalField p) (h : K.absGal) :
    ContinuousMulEquiv K.absGal K.absGal where
  toMulEquiv := MulAut.conj h
  continuous_toFun := (continuous_const.mul continuous_id).mul continuous_const
  continuous_invFun := (continuous_const.mul continuous_id).mul continuous_const

/-! ## 1. Proposition 1.1(`Interface` を取らない) -/

/-- ★**共役の経路は証明として閉じている**。

円分指標は可換群 `ℤ_[p]ˣ` への準同型なので、内部自己同型では動かない
(`χ(h⁻¹ g h) = χ(h)⁻¹ χ(g) χ(h) = χ(g)`)。
すなわち `cyclotomicCharacter_recoverable` は共役では**原理的に反証できない**
——「見つからなかった」ではなく「その経路は閉じている」ことの証明。 -/
theorem cyclotomicCharacter_conj_invariant (K : PAdicLocalField p) (h : K.absGal) :
    (cyclotomicCharacterObject (p := p)).transport (conjEquiv K h)
        ((cyclotomicCharacterObject (p := p)).obj K)
      = (cyclotomicCharacterObject (p := p)).obj K := by
  funext g
  show cyclotomicCharacter K.closure p (((MulAut.conj h).symm g).toRingEquiv)
      = cyclotomicCharacter K.closure p g.toRingEquiv
  have hs : (MulAut.conj h).symm g = h⁻¹ * g * h := by simp [MulAut.conj]
  rw [hs]
  show cyclotomicCharacter K.closure p
      ((h⁻¹ : K.absGal).toRingEquiv * g.toRingEquiv * h.toRingEquiv) = _
  rw [map_mul, map_mul]
  show (cyclotomicCharacter K.closure p (h.toRingEquiv))⁻¹ * _ * _ = _
  simp [mul_comm]

/-! ## 2. Proposition 1.2(`ResidueCardinality` を取る) -/

/-- **K′ = K の経路は自明に閉じる**。`transport` が恒等だから。

したがって `residueCard_and_degree_recoverable` の反証は
**K ≠ K′ かつ α が存在するペアの構成**にしか帰着しない。 -/
theorem residueCardAndDegree_self_route_closed (RD : ResidueCardinality p)
    (K : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K.absGal) :
    (residueCardAndDegreeObject RD).transport α ((residueCardAndDegreeObject RD).obj K)
      = (residueCardAndDegreeObject RD).obj K := rfl

/-- 我々が「2つの p進局所体が異なる」を示せる唯一の道 — 不変量(次数)の違い。 -/
theorem ne_of_finrank_ne {K K' : PAdicLocalField p}
    (h : Module.finrank ℚ_[p] K.carrier ≠ Module.finrank ℚ_[p] K'.carrier) : K ≠ K' := by
  rintro rfl; exact h rfl

/-- ★**反証は α の構成だけに帰着する**。

次数の違う2体の間に連続同型 α が1つでも作れれば、`RD` が何であっても
Proposition 1.2 の形は**即座に**反証される(`RD` の自由度すら要らない——
次数成分だけで落ちる)。

裏を返すと: **`ResidueCardinality` が自由なデータであることは、
この定理を反証する役に立たない。** 反証に要るのは α であり、
α が存在しないことこそ Proposition 1.2 が主張している内容である。
すなわち**この定理を反証するには pGC Proposition 1.2 を偽にするしかない**。 -/
theorem refutation_reduces_to_alpha (RD : ResidueCardinality p) {K K' : PAdicLocalField p}
    (hne : Module.finrank ℚ_[p] K.carrier ≠ Module.finrank ℚ_[p] K'.carrier)
    (α : ContinuousMulEquiv K.absGal K'.absGal) :
    ¬ (residueCardAndDegreeObject RD).RecoverableFromAbsGal := by
  intro h
  exact hne (congrArg Prod.snd (h α))

/-! ## 3. Corollary 1.3(`ResidueCardinality`・`SubgroupCorrespondence` を取る) -/

/-- 正規部分群は共役で動かない。 -/
theorem map_conj_of_normal {G : Type*} [Group G] (H : Subgroup G) [H.Normal] (g : G) :
    H.map (MulAut.conj g).toMonoidHom = H := by
  ext x
  simp only [Subgroup.mem_map, MulAut.conj_apply, MulEquiv.coe_toMonoidHom]
  constructor
  · rintro ⟨y, hy, rfl⟩; exact Subgroup.Normal.conj_mem ‹H.Normal› y hy g
  · intro hx
    exact ⟨g⁻¹ * x * g, by simpa using Subgroup.Normal.conj_mem ‹H.Normal› x hx g⁻¹, by group⟩

/-- ★**構成できる唯一の `SC` では、`inertia_recoverable` は証明できてしまう**。

`degenerateSC`(`Check/PGC/DischargeFiresNothing.lean`)の下では
`inertia RD SC K = ⊤` が**任意の `RD`** について成り立つ。`⊤` は任意の全射で `⊤` へ移るので、
`RecoverableFromAbsGal` が**成り立つ**。

つまり退化 witness は、Cor 3.12 のときと違って**反証を与えない**——逆に定理を与える。
`sorry` 無し。 -/
theorem inertia_recoverable_degenerateSC (RD : ResidueCardinality p) :
    (inertiaObject RD (degenerateSC (p := p))).RecoverableFromAbsGal := by
  intro K K' α
  show (inertia RD degenerateSC K).map α.toMulEquiv.toMonoidHom = inertia RD degenerateSC K'
  rw [degenerateSC_inertia_eq_top RD K, degenerateSC_inertia_eq_top RD K']
  exact Subgroup.map_top_of_surjective _ α.toMulEquiv.surjective

#print axioms cyclotomicCharacter_conj_invariant
#print axioms refutation_reduces_to_alpha
#print axioms inertia_recoverable_degenerateSC

/-!
## 探した範囲(`LeanStatus.absent (searched)` の規律)

**反証に必要なもの**は、3定理とも次のどちらかに帰着した:

- **(α) 2つの *異なる* p進局所体と、その絶対 Galois 群の間の連続同型**
- **(β) `Γ_K` の *非正規* な開部分群と、それを動かす元**

試した経路と、それぞれがどこで塞がったか:

1. **K′ = K + 恒等射** — 3件とも自明に成立。反証にならない。
2. **K′ = K + 共役**(`conjEquiv`、上で構成済み)
   - Prop 1.1: `cyclotomicCharacter_conj_invariant` で**閉じていることを証明**した。
   - Prop 1.2: `residueCardAndDegree_self_route_closed`(`transport` が恒等)。
   - Cor 1.3: `map_conj_of_normal` より、`inertia` が正規なら閉じる。
     構成できる唯一の `SC`(`degenerateSC`)では `inertia = ⊤` で正規
     ——それどころか定理そのものが証明できる(`inertia_recoverable_degenerateSC`)。
     非正規な `inertia` を作るには (β) が要る。
3. **`ULift` による同型コピー** `⟨ℚ_[p]⟩` 対 `⟨ULift ℚ_[p]⟩` — **2段階で塞がった**:
   - mathlib に **`Field (ULift α)` インスタンスが無い**(2026-08-14 実測、
     `inferInstance` が失敗)。よって `PAdicLocalField p` の値として構成できない。
   - 仮にインスタンスを足しても、`ℚ_[p] ≠ ULift ℚ_[p]`(**型の非等号**)は
     Lean で証明できない。同型な型を区別する手段が無い。
4. **次数の違う2体**(例: `ℚ_[p]` と2次拡大)— こちらは `ne_of_finrank_ne` で
   **異なることは示せる**が、その間の α は存在しない
   (存在すれば `refutation_reduces_to_alpha` で Prop 1.2 が落ちる)。
   すなわち**この経路は原理的に反証に使えない**。
5. **非正規な開部分群の直接構成** — `Γ_{ℚ_p}` の部分群構造の知識が要る。
   指数 2 の部分群は常に正規なので、非 Galois な3次拡大以上が要る。
   `AlgebraicClosure ℚ_[p]` の中にそれを構成する道具は無い(mathlib 実測: p進体の
   代数拡大の分類・分岐理論は 0 件——PLAN 事実3)。

**まとめ**: 3定理の反証は、いずれも「我々がまだ作れない数学的対象」に帰着した。
これは `residueCard base = p` を越えて f ≥ 2 の実例が作れなかったのと**同じ壁**である
(`Check/IUTchIII/...` ではなく `Check/PGC/DischargeFiresNothing.lean` の未解決節を参照)。

## 判断への含意

Cor 3.12 が反証できて pGC の3定理が反証できないのは、定理の質の差ではなく
**`Interface` が抽象的な型の上にあるか、具体的な mathlib 対象の上にあるか**の差である。
IUT 方向の `Interface` は(mathlib に望月語彙が 0 件である以上)前者にならざるを得ない。
-/

end ABC3.Check.PGC
