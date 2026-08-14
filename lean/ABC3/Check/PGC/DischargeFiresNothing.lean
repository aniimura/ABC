import ABC3.Check.PGC.InertiaDegeneracyMoved
import ABC3.Found.PGC.ResidueCardinality

/-!
# discharge しても検査は発火しなかった — PLAN §3 の主張の実測

PLAN §3 はこう主張している:

> Track B が `Interface` の実例を1つ供給した瞬間、それに依存する Track A の
> statement 群が一斉に非空虚性の検査を受ける。**空虚さが後から自動的に暴かれる構造**
> になっている

2026-08-14 に `ResidueCardinality` を discharge した(`Found/PGC/ResidueCardinality.lean`)。
**では実際に何か暴かれたか。** このファイルはその問いを機械にかける。

答えは **暴かれなかった**。理由は2つあり、どちらも下で証明する。

## (i) 退化は `RD` から `SC` へ移っただけだった

`Check/PGC/InertiaDegeneracyMoved.lean` は「偽の `RD` を入れると `inertia = ⊤`」を示した。
本物の `RD` を入れれば排除される、というのが期待だった。**排除されない。**

`degenerateSC`(どの開部分群にも `K` 自身を返す `SubgroupCorrespondence`)を取ると、
`degenerateSC_inertia_eq_top` が示すとおり **`RD` が何であっても** `inertia = ⊤` に潰れる。
本物の `realResidueCardinality` も例外ではない(`realRD_inertia_eq_top`)。

しかも潰れる理由は `Interface` の条件そのものにある——`isPrimePow` が
`RD.card K = p ^ f`(`f > 0`)を要求するので `RD.card K > 1` が常に成り立ち、
`q = q ^ [Γ_K : H]` は `[Γ_K : H] = 1` を強制する。
**`RD` を本物にしても、この経路は塞げない。** 塞ぐ情報は `SC` の側にある。

## (ii) G2 の `Nonempty` は退化 witness で満たせる

そもそも `Nonempty (ResidueCardinality p)` は Track B を待たずに示せた——
`degenerateRD` で足りる(`residueCardinality_nonempty_by_degenerate`)。
つまり `ResidueCardinality` を仮説に持つ statement 群は、
**discharge 前から空虚ではなかった**。G2 が検査しているのは「実例があるか」であって
「その実例が意図したものか」ではない。

同じことが、まだ待ち行列に残っている `SubgroupCorrespondence` にも言える
(`subgroupCorrespondence_nonempty_by_degenerate`)。
★そちらは**意図的に `SubgroupCorrespondence.nonvacuous` と名付けていない**——
check.mjs の G2 は宣言名で witness を探すので、その名前を付けると
退化 witness で作業キューが空になる。本プロジェクトが最も警戒している失敗の形。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## (i) 退化は `SC` へ移った -/

/-- `Interface` の条件だけから、剰余体の位数は必ず 1 より大きい。

`isPrimePow` の `0 < f` がここで効いている。 -/
theorem one_lt_card (RD : ResidueCardinality p) (K : PAdicLocalField p) : 1 < RD.card K := by
  obtain ⟨f, hf, hc⟩ := RD.isPrimePow K
  rw [hc]
  exact Nat.one_lt_pow hf.ne' (Fact.out : p.Prime).one_lt

/-- 退化した `SubgroupCorrespondence`: どの開部分群にも `K` 自身を返す。

`field_top`(H = ⊤ なら L = K)は `rfl` で満たされる——
`Interface` が `SC` に課している条件は**これだけ**なので、この定義は合法である。 -/
def degenerateSC : SubgroupCorrespondence p where
  field := fun K _ _ => K
  field_top := fun _ _ => rfl

/-- ★**本題**: `RD` が何であっても、退化した `SC` の下では惰性群は `⊤` に潰れる。

`degenerateRD_inertia_eq_top`(`InertiaDegeneracyMoved.lean`)は
「偽の `RD` × 任意の `SC`」で潰れることを示していた。こちらは
「**任意の `RD`** × 偽の `SC`」で潰れることを示す。退化への経路は2本あり、
`RD` を本物にしても塞がるのは片方だけ——実際には片方も塞がっていない(下記)。 -/
theorem degenerateSC_inertia_eq_top (RD : ResidueCardinality p) (K : PAdicLocalField p) :
    inertia RD degenerateSC K = ⊤ := by
  have hsub : {H : Subgroup K.absGal |
      ∃ hH : IsOpen (H : Set K.absGal), IsUnramifiedAt RD degenerateSC K H hH}
      ⊆ {⊤} := by
    rintro H ⟨hH, hEq⟩
    -- `degenerateSC.field K H hH = K` なので、判定条件は `q = q ^ [Γ_K : H]` に潰れる
    have h1 : RD.card K = RD.card K ^ H.index := by
      simpa [degenerateSC, IsUnramifiedAt] using hEq
    have h1' : RD.card K ^ 1 = RD.card K ^ H.index := by rw [pow_one]; exact h1
    have h2 : (1 : ℕ) = H.index := Nat.pow_right_injective (one_lt_card RD K) h1'
    exact Set.mem_singleton_iff.mpr (Subgroup.index_eq_one.mp h2.symm)
  have h4 : (⊤ : Subgroup K.absGal) ≤ inertia RD degenerateSC K := by
    have h5 := sInf_le_sInf hsub
    rwa [sInf_singleton] at h5
  exact top_le_iff.mp h4

/-- ★**問1への答え**: 本物の `realResidueCardinality` を入れても退化は排除されない。 -/
theorem realRD_inertia_eq_top (K : PAdicLocalField p) :
    inertia (realResidueCardinality p) degenerateSC K = ⊤ :=
  degenerateSC_inertia_eq_top _ K

/-- 対照: 偽の `RD` でも同じ結論。**両者は `degenerateSC` の上で区別できない。** -/
theorem degenerateRD_inertia_eq_top' (K : PAdicLocalField p) :
    inertia (degenerateRD (p := p)) degenerateSC K = ⊤ :=
  degenerateSC_inertia_eq_top _ K

/-! ## (i-b) そもそも本物と偽物は区別できるのか

問2: `realResidueCardinality` と `degenerateRD` は**関数として異なるか**。

`ResidueCardinality` の第2フィールドは `Prop` なので、2つの値が等しいことは
`card` が等しいことと同値(`eq_iff_card_eq`)。したがって問は
「剰余次数 f > 1 の p進局所体が存在するか」に**厳密に還元される**
(`real_ne_degenerate_iff`)。

数学的には存在する(不分岐拡大)。しかし Lean 上では、現在
`PAdicLocalField p` の実例として名指しできるのは `base`(= `ℚ_[p]` 自身、f = 1)だけで、
この同値の右辺を示す手段が無い。詳細は本ファイル末尾の「未解決」節。 -/

/-- `ResidueCardinality` は `card` フィールドだけで決まる(`isPrimePow` は `Prop`)。 -/
theorem eq_iff_card_eq {r₁ r₂ : ResidueCardinality p} : r₁ = r₂ ↔ r₁.card = r₂.card := by
  constructor
  · rintro rfl; rfl
  · cases r₁; cases r₂; rintro h; cases h; rfl

/-- ★**問2の還元**: 本物と偽物が異なることは、
剰余次数 f > 1 の p進局所体が存在することと**同値**。 -/
theorem real_ne_degenerate_iff :
    realResidueCardinality p ≠ degenerateRD ↔ ∃ K : PAdicLocalField p, residueCard K ≠ p := by
  rw [Ne, eq_iff_card_eq]
  show ¬ (residueCard = fun _ : PAdicLocalField p => p) ↔ _
  rw [funext_iff]
  push Not
  rfl

/-! ## (ii) G2 の `Nonempty` は退化 witness で満たせる -/

/-- ★`Nonempty (ResidueCardinality p)` は Track B の実装を待たずに示せた。

したがって「Track B が実例を供給した瞬間に空虚さが暴かれる」は、
少なくとも `ResidueCardinality` については**起こりようがなかった**——
実例は最初から(`degenerateRD` として)存在していた。 -/
theorem residueCardinality_nonempty_by_degenerate : Nonempty (ResidueCardinality p) :=
  ⟨degenerateRD⟩

/-- 同じことが、まだ待ち行列に残る `SubgroupCorrespondence` にも言える。

★**この宣言を `SubgroupCorrespondence.nonvacuous` と名付けてはならない**——
check.mjs の G2 は宣言名で witness を探すので、その名前を付けた瞬間に
`SubgroupCorrespondence` は作業キューから消える。実装は何も進んでいないのに。 -/
theorem subgroupCorrespondence_nonempty_by_degenerate : Nonempty (SubgroupCorrespondence p) :=
  ⟨degenerateSC⟩

/-!
## 未解決 — `real_ne_degenerate_iff` の右辺は示せていない

`∃ K : PAdicLocalField p, residueCard K ≠ p` は**数学的には真**(ℚ_p は各次数の
不分岐拡大を持つ)。しかし Lean 上では示せていない。何が要るかを実測して記録する。

### 何が足りないか(2026-08-14 実測)

必要なのは剰余次数 f ≥ 2 の p進局所体だが、それには2段が要る:

1. **不分岐拡大そのもの** — mathlib v4.31.0-rc2 に構成が無い。
   `Mathlib/NumberTheory/LocalField/` は `Basic.lean` 1ファイルのみ。
   `grep -rliE "unramified" Mathlib/NumberTheory/ Mathlib/RingTheory/Valuation/` が返すのは
   `NumberField/*`(大域)と `RamificationInertia/*`(Dedekind 環の一般論)で、
   「完備離散付値体は各次数の不分岐拡大を持つ」という形の結果は無い。
   `grep -rn "unramified extension"` → `NumberField/ExistsRamified.lean` の1件のみ(別の主張)。
   なお**拡大体を作るだけ**なら `AdjoinRoot` で足りるが、それだけでは `residueCard` を
   計算できないので用を成さない。剰余次数を出すには不分岐性の理論そのものが要る。

2. **我々の `residueCard` に繋ぐ橋** — 上が仮にあっても、`residueCard` は
   `Found/PGC/LocalFieldNorm.lean` の**スペクトルノルム由来の** `Valued` 構造に対する
   `Nat.card 𝓀[·]` なので、`Ideal.inertiaDeg` 等の代数的な剰余次数と繋ぐ必要がある。

### 一番簡単な場合(`base` = ℚ_p、f = 1)ですら未完了

`residueCard base = p`(= 我々の唯一の実例の上で本物と偽物が一致すること)を
証明しようとして、2つの具体的な障害に当たった(2026-08-14、時間を区切って中断):

- **障害A(ノルムのダイヤモンド、再来)**: `ℚ_[p]` 上に `NormedField` が2つある——
  標準のものと、`(selfField p).carrier` に付くスペクトルノルム由来のもの。
  `x : (selfField p).carrier` に対する `‖x‖` は後者に解決されるため、
  `PadicInt.mem_subring_iff`(`x ∈ subring p ↔ ‖x‖ ≤ 1`、標準ノルム)と
  **定義的に一致しない**。型注釈 `(x : ℚ_[p])` ではインスタンス選択は変わらない。
  両者が値として等しいことは `norm_algebraMap` から出るが、`show`/`rw` の
  段階で syntactic に噛み合わせる作業が要る。
- **障害B**: `IsLocalRing ↥(PadicInt.subring p)` がインスタンス探索で見つからない。
  mathlib が持つのは `IsLocalRing ℤ_[p]` で、`↥(PadicInt.subring p)` とは
  定義的に等しいが構文的に別。`PadicInt.residueField`
  (`IsLocalRing.ResidueField ℤ_[p] ≃+* ZMod p`)へ繋ぐには明示的な転送が要る。

材料は揃っている(`PadicInt.subring` / `PadicInt.residueField` /
`IsLocalRing.ResidueField.mapEquiv` / `RingEquiv.subringCongr`)。
足りないのは上記2つの噛み合わせだけで、**数学ではなく配線の問題**である。
-/

end ABC3.Check.PGC
