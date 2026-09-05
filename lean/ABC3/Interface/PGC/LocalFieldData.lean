import ABC3.Skeleton.PGC.Setup

/-!
# [pGC] §1 — まだ構築できていない基礎

原典が §1 で**所与として**使うが、我々がまだ mathlib から得られないもの。

`axiom` ではなく `structure` で受ける(`ABC3/Meta/Calibration.lean` の実演を参照)。
各 `structure` は非空虚 witness を持つか、何を待っているかを書く(G2)。

## ★import の向き — `Interface` は `Found` を import しない(2026-08-14 に確定)

一度、G2 の witness を `structure` と同じファイルに置くために
本ファイルから `Found/PGC/LocalFieldNorm` を import した。**これは誤りだった**。

`Skeleton/PGC/Section1.lean` は本ファイルを import するので、その向きだと
**`Skeleton` が `Found` を推移的に import する**ことになる。`Skeleton` が
`(RD : ResidueCardinality p)` を仮説に取るのは「実装が無くても statement を書ける」
ためであり、実装に依存しないことが要点(PLAN §3 の2トラック構成そのもの)。
その向きでは実装が無いと `Skeleton` がビルドできず、条件付き形式化の設計が壊れる。

**規則**: `Interface/` から `Found/` を import してはならない。
`check.mjs` が検査する(fixture D23/D24)。

したがって非空虚 witness は実装側に置く——
`ResidueCardinality.nonvacuous` は `Found/PGC/ResidueCardinality.lean` にある。
check.mjs の G2 は宣言名を木全体から探すので、この配置で通る。
-/

namespace ABC3.Interface.PGC

open ABC3.Meta ABC3.Skeleton.PGC

variable (p : ℕ) [Fact p.Prime]

/-- 各 p進局所体 K に、剰余体の元の個数 q を与えるデータ。

原文 [pGC] 物理 p.3:

原文 (pGC p.3):
> Let k be the residue field of O[scr]_K (the ring of integers of K). Thus, k is the
> field of q = p^f elements.

**なぜ Interface なのか(実測、2026-08-14)**:
mathlib は `IsNonarchimedeanLocalField`(剰余体の有限性込み)を持つが、
**「`ℚ_[p]` の有限次拡大」からそのインスタンスを導く経路が無い**——
`IsNonarchimedeanLocalField` は `[ValuativeRel K] [TopologicalSpace K]` を
入力として要求し、有限次拡大に付値構造を与える部分が mathlib に無い。

公開プロジェクト `kbuzzard/ClassFieldTheory` の `IsNonarchimedeanLocalField/Instances.lean`
が同じ穴を埋めようとしているが、**`sorry` が 11 件残る**(測定日 2026-08-14、
記録は `ResearchPaper/lean-ecosystem.json`)。

**★2026-08-14: discharge 済み**。上の「なぜ Interface なのか」は
`IsNonarchimedeanLocalField` へ繋ぐ経路についての測定としては今も正しいが、
**この `structure` を埋めるのにその経路は要らなかった**——`ℚ_[p]` 上の
スペクトルノルムで直接 `NormedField`/`Valued` を入れれば、剰余体の有限性
(`Found/ResidueFieldFinite.lean`)と標数(`CharP.charP_iff_prime_eq_zero`)が
そのまま出る。`waiting` に書いていた
「`IsNonarchimedeanLocalField` へ繋ぐ」は**必要条件ではなかった**。
実装は `Found/PGC/LocalFieldNorm.lean` と `Found/PGC/ResidueCardinality.lean`。

**G2 の非空虚 witness `ResidueCardinality.nonvacuous` は
`Found/PGC/ResidueCardinality.lean` にある**(上記「import の向き」参照)。
ここに置くと `Skeleton` が `Found` を推移的に import してしまう。 -/
structure ResidueCardinality where
  /-- 剰余体の元の個数 q -/
  card : PAdicLocalField p → ℕ
  /-- 原文「k is the field of q = p^f elements」——q は p の正の冪。
      この条件が無いと `card := fun _ => 0` でも通ってしまい、内容が消える。
      この主張自体の検査は `Check/PGC/ResidueCardinalityNondegenerate.lean`。 -/
  isPrimePow : ∀ K, ∃ f : ℕ, 0 < f ∧ card K = p ^ f

/-- 開部分群 H ⊆ Γ_K に対応する中間体 L もまた p進局所体である、という対応。

原文 (pGC p.3):
> Now let H ⊆ Γ_K be an open subgroup. Let L ⊇ K be the extension field of K
> corresponding to H. By applying Proposition 1.2 to L and H, we see that the number
> q_L of elements in the residue field of O[scr]_L can be recovered group-theoretically
> from H ⊆ Γ_K.

原文は「Proposition 1.2 を (L, H) に適用する」と述べており、
**L が p進局所体であり、その絶対 Galois 群が H である**ことを前提している。
その前提をここで型にする。

**なぜ Interface なのか**: mathlib は Galois 対応(`IntermediateField.fixedField`)と
Krull 位相を持つが、「開部分群に対応する中間体が ℚ_p 上有限次である」ことと
「その絶対 Galois 群が H と同一視できる」ことは、この設定では未接続。

★`field_top` は原文が明示していない——H = Γ_K のとき L = K であることは
対応の定義から従う自明な事実だが、**我々が明示的に課している**ことを記録しておく。
これが無いと不分岐な開部分群の集合が空になりうる。

**★★★2026-09-05: discharge 済み**。`waiting` が挙げていた二つ

1. 開部分群 H に対応する中間体が p進局所体であること
2. その絶対 Galois 群が H であること

は両方とも構成された:
1. `Found/PGC/SubgroupCorrespondenceConstruction.lean::fixedFieldLocalField`
   (Krull 位相の開部分群の固定体は ℚ_p 上有限次——
   `krullTopology_mem_nhds_one_iff` で有限次中間体に挟む)
2. `Found/PGC/AdjoinFieldClosure.lean::absGalFixedFieldEquiv`
   (代数閉包の同一視 + `IntermediateField.fixingSubgroupEquiv` +
   `InfiniteGalois.fixingSubgroup_fixedField`)

実物は `Found/PGC/SubgroupCorrespondenceConstruction.lean::subgroupCorrespondence`。
非空虚 witness は同ファイルの `SubgroupCorrespondence.nonvacuous`
(`Interface` は `Found` を import できないので、実装側から名前空間に足す)。 -/
structure SubgroupCorrespondence where
  /-- 開部分群 H に対応する中間体(これも p進局所体) -/
  field : (K : PAdicLocalField p) → (H : Subgroup K.absGal) →
    IsOpen (H : Set K.absGal) → PAdicLocalField p
  /-- H = Γ_K に対応するのは K 自身 -/
  field_top : ∀ K h, field K ⊤ h = K

def SubgroupCorrespondence.waiting : WaitingFor :=
  { what := "開部分群 H ⊆ Γ_K に対応する中間体が p進局所体であり、その絶対 Galois 群が H であること"
    trackB := "Found/LocalField — Krull 位相の Galois 対応と、ℚ_p 上の有限次性の接続" }

/-- **[pGC] §2** Proposition 2.2 以降が仮説として使う、高次分岐群(上付き番号付け)の族
`{Γ_K^v}_{v>0}`。

原文 (pGC p.4):
> Then we shall denote by Γ^v_K ⊆ Γ_K the higher ramification group associated to the
> number v in the "upper numbering" (see, e.g., [3], p. 155).

**なぜ Interface なのか**: これを構成するには Herbrand の定理
(Theorem 1 of [3], p.155: 上付き番号付けの分岐群の Γ_K^ab での像が U_K^v に一致する)が
要る。mathlib に**不在**(実測: `Mathlib/RingTheory/Valuation/RamificationGroup.lean` に
`TODO: Define higher ramification groups in lower numbering` と明記、上付き番号付けは
それ以前の段階)。

`Gv` の定義域は原文では `v > 0` のみだが、`Antitone` の条件を書きやすくするため
全実数上で定義する(`v ≤ 0` での値を固定する、という我々の設計判断——逸脱として記録)。
`Found/PGC/FilteredGroup.lean` の `FilteredGroup` の3条件(閉・正規・反単調)と同形。 -/
structure RamificationFiltration (p : ℕ) [Fact p.Prime] where
  /-- `Γ_K^v ⊆ Γ_K`(v は実数) -/
  Gv : (K : PAdicLocalField p) → ℝ → Subgroup K.absGal
  /-- 各 `Γ_K^v` は閉部分群 -/
  isClosed : ∀ K v, IsClosed (Gv K v : Set K.absGal)
  /-- 各 `Γ_K^v` は正規部分群 -/
  isNormal : ∀ K v, (Gv K v).Normal
  /-- 降下条件: `v1 ≥ v2 → Γ_K^{v1} ⊆ Γ_K^{v2}` -/
  antitone : ∀ K, Antitone (Gv K)

def RamificationFiltration.waiting (p : ℕ) [Fact p.Prime] : WaitingFor :=
  { what := "高次分岐群(上付き番号付け){Γ_K^v}_{v>0} の族。Herbrand の定理を要する"
    trackB := "Found/PGC — RamificationGroup.lean の上付き/下付き番号付けの拡張(mathlib 本体の TODO)" }

end ABC3.Interface.PGC
