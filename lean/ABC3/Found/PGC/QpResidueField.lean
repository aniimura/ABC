import ABC3.Found.PGC.LocalFieldNorm
import Mathlib.NumberTheory.Padics.RingHoms

/-!
# `ℚ_[p]` の剰余体は `p` 元体 — `residueCard` の**値**を初めて計算する

`Found/PGC/LocalFieldNorm.lean` は `residueCard K = Nat.card 𝓀[K]` を定義し、
それが `p` の正の冪であることまで示した。しかし **値を1つも計算していなかった**。
計算できなければ、`realResidueCardinality` が本物の数学に繋がっている保証は無く、
退化した `degenerateRD` に対する優位も我々の側からは示せない
(`Check/PGC/DischargeFiresNothing.lean` の議論)。

ここで最初の値を出す: **`residueCard (selfField p) = p`**(`selfField p` = `ℚ_[p]` 自身)。

## ★ダイヤモンドの避け方(2026-08-14。前回はここで止まった)

`ℚ_[p]` の上には `NormedField` が **2つ** ある——標準のものと、
`(selfField p).carrier` に付くスペクトルノルム由来のもの。前回はこの2つを
噛み合わせようとして止まった。今回は**噛み合わせずに済ませた**。

鍵は3つ:

1. **ノルムの一致は `algebraMap` 経由で述べる。**
   `‖y‖ = ‖y‖` と書くとどちらのインスタンスに解決されるか制御できない
   (型注釈 `(y : ℚ_[p])` では変わらない——項の型は変わらないから)。
   既存の `norm_algebraMap K x : ‖algebraMap ℚ_[p] K.carrier x‖ = ‖x‖` は
   **左右で型が違う**ので、インスタンスが構文的に確定している。
   `selfField` に適用すると、そのまま欲しい等式になる(`norm_selfField`)。
2. **ノルム構造の一致ではなく、付値環の集合としての一致を示す。**
   `integer_selfField : 𝒪[(selfField p).carrier] = PadicInt.subring p`(`Subring ℚ_[p]` の等式)。
   ここまで来れば、以降にノルムは出てこない。
3. **中間型 `↥(PadicInt.subring p)` を型として書かない。**
   前回 `IsLocalRing ↥(PadicInt.subring p)` が見つからず止まったが、
   `RingEquiv.subringCongr` の**行き先を `ℤ_[p]` と宣言してしまえば**
   (両者は定義的に等しい)、インスタンス探索は素直に `IsLocalRing ℤ_[p]` を見つける。
   `letI` で局所インスタンスを与える方法は失敗する——`ResidueField` が
   インスタンスに依存するので、局所 fvar 版と `PadicInt` 版が単一化しない。

`Found/ResidueFieldFinite.lean` の教訓「ダイヤモンドは**触れずに済ませる**方が安い」の
2例目。前回の中断報告では「障害A・障害B」と書いたが、**どちらも正面から解く必要は無かった**。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable (p : ℕ) [Fact p.Prime]

/-- `ℚ_[p]` 自身を p進局所体と見たもの(有限次拡大の次数 1 の場合)。 -/
noncomputable def selfField : PAdicLocalField p := { carrier := ℚ_[p] }

/-- スペクトルノルムは `ℚ_[p]` 上で標準の p進ノルムに一致する。

★左辺のノルムは `(selfField p).carrier` 上のスペクトルノルム、
右辺は `ℚ_[p]` 上の標準ノルム。**型が違うのでインスタンスが混ざらない**。 -/
theorem norm_selfField (y : (selfField p).carrier) : ‖y‖ = @norm ℚ_[p] _ y :=
  norm_algebraMap (selfField p) y

open scoped NormedField Valued in
/-- 付値環は `ℤ_[p]`(を `ℚ_[p]` の部分環と見たもの)にほかならない。

ここから先にノルムは現れない——**ダイヤモンドはここで置き去りにする**。 -/
theorem integer_selfField :
    (𝒪[(selfField p).carrier] : Subring ℚ_[p]) = PadicInt.subring p := by
  ext y
  rw [Valued.integer.mem_iff]
  show (Valued.v y ≤ 1) ↔ (@norm ℚ_[p] _ y ≤ 1)
  rw [← norm_selfField p y]
  simp [← NNReal.coe_le_coe]

open scoped NormedField Valued in
/-- 付値環と `ℤ_[p]` の環同型。

★行き先を `↥(PadicInt.subring p)` ではなく **`ℤ_[p]` と宣言する**のが要点
(両者は定義的に等しい。mathlib の `PadicInt.instCommRing` 自身が
`inferInstanceAs <| CommRing (subring p)` で定義されている)。 -/
noncomputable def integerEquivPadicInt : ↥(𝒪[(selfField p).carrier]) ≃+* ℤ_[p] :=
  RingEquiv.subringCongr (integer_selfField p)

open scoped NormedField Valued in
/-- 剰余体は `ZMod p`。 -/
noncomputable def residueFieldEquivZMod : 𝓀[(selfField p).carrier] ≃+* ZMod p :=
  (IsLocalRing.ResidueField.mapEquiv (integerEquivPadicInt p)).trans PadicInt.residueField

/-- ★**`residueCard` の値を計算した最初の例**: `ℚ_[p]` の剰余体は `p` 元体。

これで `realResidueCardinality` は「型が付く」だけの対象ではなくなった——
少なくとも1点で値が出る。 -/
theorem residueCard_selfField : residueCard (selfField p) = p := by
  rw [residueCard, Nat.card_congr (residueFieldEquivZMod p).toEquiv, Nat.card_zmod]

end ABC3.Found.PGC
