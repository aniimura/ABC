import ABC3.Skeleton.PGC.Setup

/-!
# [pGC] §1 — まだ構築できていない基礎

原典が §1 で**所与として**使うが、我々がまだ mathlib から得られないもの。

`axiom` ではなく `structure` で受ける(`ABC3/Meta/Calibration.lean` の実演を参照)。
各 `structure` は非空虚 witness を持つか、何を待っているかを書く(G2)。
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

**非空虚性について**: この `structure` は実際には充足可能である(K は局所体なので
剰余体は有限で位数は p の冪)。ただし**その事実を我々はまだ証明できない**ので、
非空虚 witness ではなく `waiting` を置く。 -/
structure ResidueCardinality where
  /-- 剰余体の元の個数 q -/
  card : PAdicLocalField p → ℕ
  /-- 原文「k is the field of q = p^f elements」——q は p の正の冪。
      この条件が無いと `card := fun _ => 0` でも通ってしまい、内容が消える。 -/
  isPrimePow : ∀ K, ∃ f : ℕ, 0 < f ∧ card K = p ^ f

def ResidueCardinality.waiting : WaitingFor :=
  { what := "ℚ_p の有限次拡大 K の整数環 𝒪_K と、その剰余体が位数 p^f の有限体であること"
    trackB := "Found/LocalField — 有限次拡大に付値構造を与え IsNonarchimedeanLocalField へ繋ぐ" }

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
これが無いと不分岐な開部分群の集合が空になりうる。 -/
structure SubgroupCorrespondence where
  /-- 開部分群 H に対応する中間体(これも p進局所体) -/
  field : (K : PAdicLocalField p) → (H : Subgroup K.absGal) →
    IsOpen (H : Set K.absGal) → PAdicLocalField p
  /-- H = Γ_K に対応するのは K 自身 -/
  field_top : ∀ K h, field K ⊤ h = K

def SubgroupCorrespondence.waiting : WaitingFor :=
  { what := "開部分群 H ⊆ Γ_K に対応する中間体が p進局所体であり、その絶対 Galois 群が H であること"
    trackB := "Found/LocalField — Krull 位相の Galois 対応と、ℚ_p 上の有限次性の接続" }

end ABC3.Interface.PGC
